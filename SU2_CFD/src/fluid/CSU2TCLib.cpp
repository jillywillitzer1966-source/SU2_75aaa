/*!
 * \file CSU2TCLib.cpp
 * \brief Source of user defined 2T nonequilibrium gas model.
 * \author C. Garbacz, W. Maier, S. R. Copeland, J. Needels
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/fluid/CSU2TCLib.hpp"
#include "../../../Common/include/option_structure.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <unordered_map>

namespace {

struct CustomReactionRow {
  unsigned short reaction_index = 0;
  unsigned short reactants[3] = {0, 0, 0};
  unsigned short products[3] = {0, 0, 0};
  su2double arrhenius_a = 0.0;
  su2double arrhenius_eta = 0.0;
  su2double arrhenius_theta = 0.0;
  su2double tcf_a = 1.0;
  su2double tcf_b = 0.0;
  su2double tcb_a = 1.0;
  su2double tcb_b = 0.0;
  unsigned short keq_table_id = 0;
};

static string Trim(const string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == string::npos) return "";
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

static vector<string> SplitCsvLikeLine(string line) {
  for (char& c : line) {
    if (c == '\t') c = ',';
  }

  vector<string> fields;
  string current;
  bool in_quotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        current.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
      continue;
    }
    if (c == ',' && !in_quotes) {
      fields.push_back(Trim(current));
      current.clear();
      continue;
    }
    current.push_back(c);
  }
  fields.push_back(Trim(current));
  return fields;
}

static string CanonicalHeader(string s) {
  string out;
  out.reserve(s.size());
  for (unsigned char ch : s) {
    if (std::isalnum(ch)) out.push_back(static_cast<char>(std::tolower(ch)));
    else out.push_back('_');
  }
  return out;
}

static int FindHeaderByTokens(const vector<string>& header, const vector<string>& tokens) {
  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    const string canon = CanonicalHeader(header[i]);
    bool match = true;
    for (const auto& token : tokens) {
      if (canon.find(token) == string::npos) {
        match = false;
        break;
      }
    }
    if (match) return i;
  }
  return -1;
}

static bool TryParseInt(const string& s, int& value) {
  const string t = Trim(s);
  if (t.empty()) return false;
  char* end_ptr = nullptr;
  const long parsed = std::strtol(t.c_str(), &end_ptr, 10);
  if (end_ptr == t.c_str() || *end_ptr != '\0') return false;
  value = static_cast<int>(parsed);
  return true;
}

static bool TryParseDouble(const string& s, su2double& value) {
  const string t = Trim(s);
  if (t.empty()) return false;
  char* end_ptr = nullptr;
  const double parsed = std::strtod(t.c_str(), &end_ptr);
  if (end_ptr == t.c_str() || *end_ptr != '\0') return false;
  value = parsed;
  return true;
}

static unsigned short ParseSpeciesSlot(const string& s, unsigned short nSpecies, const string& context) {
  int value = -1;
  if (!TryParseInt(s, value)) return nSpecies;
  if (value == static_cast<int>(nSpecies)) return nSpecies;
  if (value < 0 || value >= static_cast<int>(nSpecies)) {
    SU2_MPI::Error("CONFIG ERROR: Species index out of range in " + context, CURRENT_FUNCTION);
  }
  return static_cast<unsigned short>(value);
}

} // namespace

void CSU2TCLib::LoadCustomChemistryTables(const CConfig* config, bool viscous) {

  const string species_file = config->GetSU2_NONEQ_SpeciesTable_File();
  const string electronic_file = config->GetSU2_NONEQ_ElectronicStates_File();
  const string reactions_file = config->GetSU2_NONEQ_Reactions_File();
  const string transport_file = config->GetSU2_NONEQ_Transport_File();
  const string keq_file = config->GetSU2_NONEQ_KeqCoeff_File();

  auto read_table = [](const string& file_name) {
    ifstream in(file_name.c_str());
    if (!in) {
      SU2_MPI::Error("CONFIG ERROR: Cannot open custom chemistry table file: " + file_name, CURRENT_FUNCTION);
    }

    vector<vector<string>> rows;
    string line;
    while (std::getline(in, line)) {
      if (Trim(line).empty()) continue;
      rows.push_back(SplitCsvLikeLine(line));
    }

    if (rows.empty()) {
      SU2_MPI::Error("CONFIG ERROR: Empty custom chemistry table file: " + file_name, CURRENT_FUNCTION);
    }
    return rows;
  };

  /*--- Species table ---*/
  {
    const auto rows = read_table(species_file);
    const auto& header = rows[0];

    const int c_idx = FindHeaderByTokens(header, {"index"});
    const int c_mw = FindHeaderByTokens(header, {"molar", "mass"});
    const int c_rot = FindHeaderByTokens(header, {"rotation", "modes"});
    const int c_hf = FindHeaderByTokens(header, {"formation", "enthalpy"});
    const int c_ref = FindHeaderByTokens(header, {"reference", "temperature"});
    const int c_tv = FindHeaderByTokens(header, {"char", "vib", "temperature"});

    if (c_idx < 0 || c_mw < 0 || c_rot < 0 || c_hf < 0 || c_ref < 0 || c_tv < 0) {
      SU2_MPI::Error("CONFIG ERROR: Missing required columns in species table.", CURRENT_FUNCTION);
    }

    vector<bool> seen(nSpecies, false);

    for (size_t iRow = 1; iRow < rows.size(); ++iRow) {
      const auto& row = rows[iRow];
      if (c_idx >= static_cast<int>(row.size())) continue;

      int idx = -1;
      if (!TryParseInt(row[c_idx], idx)) continue;
      if (idx < 0 || idx >= static_cast<int>(nSpecies)) {
        SU2_MPI::Error("CONFIG ERROR: Species index out of range in species table.", CURRENT_FUNCTION);
      }

      su2double mw = 0.0, rot = 0.0, hf = 0.0, tref = 0.0, tv = 0.0;
      if (c_mw >= static_cast<int>(row.size()) || !TryParseDouble(row[c_mw], mw)) {
        SU2_MPI::Error("CONFIG ERROR: Missing molar mass value in species table.", CURRENT_FUNCTION);
      }
      if (c_rot >= static_cast<int>(row.size()) || !TryParseDouble(row[c_rot], rot)) {
        SU2_MPI::Error("CONFIG ERROR: Missing rotation modes value in species table.", CURRENT_FUNCTION);
      }
      if (c_hf >= static_cast<int>(row.size()) || !TryParseDouble(row[c_hf], hf)) {
        SU2_MPI::Error("CONFIG ERROR: Missing formation enthalpy value in species table.", CURRENT_FUNCTION);
      }
      if (c_ref >= static_cast<int>(row.size()) || !TryParseDouble(row[c_ref], tref)) {
        SU2_MPI::Error("CONFIG ERROR: Missing reference temperature value in species table.", CURRENT_FUNCTION);
      }
      if (c_tv >= static_cast<int>(row.size()) || !TryParseDouble(row[c_tv], tv)) {
        SU2_MPI::Error("CONFIG ERROR: Missing vibrational characteristic temperature value in species table.", CURRENT_FUNCTION);
      }

      MolarMass[idx] = mw;
      RotationModes[idx] = rot;
      Enthalpy_Formation[idx] = hf;
      Ref_Temperature[idx] = tref;
      CharVibTemp[idx] = tv;
      seen[idx] = true;
    }

    for (unsigned short iSp = 0; iSp < nSpecies; ++iSp) {
      if (!seen[iSp]) {
        SU2_MPI::Error("CONFIG ERROR: Species table does not provide contiguous species entries from 0 to nSpecies-1.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Electronic states table ---*/
  {
    const auto rows = read_table(electronic_file);
    const auto& header = rows[0];

    const int c_sp = FindHeaderByTokens(header, {"species", "index"});
    const int c_lv = FindHeaderByTokens(header, {"level", "index"});
    const int c_te = FindHeaderByTokens(header, {"char", "el", "temperature"});
    const int c_g = FindHeaderByTokens(header, {"degeneracy"});

    if (c_sp < 0 || c_lv < 0 || c_te < 0 || c_g < 0) {
      SU2_MPI::Error("CONFIG ERROR: Missing required columns in electronic states table.", CURRENT_FUNCTION);
    }

    vector<map<unsigned short, pair<su2double, su2double>>> state_map(nSpecies);

    for (size_t iRow = 1; iRow < rows.size(); ++iRow) {
      const auto& row = rows[iRow];
      if (c_sp >= static_cast<int>(row.size()) || c_lv >= static_cast<int>(row.size())) continue;

      int iSp = -1, iLv = -1;
      if (!TryParseInt(row[c_sp], iSp) || !TryParseInt(row[c_lv], iLv)) continue;
      if (iSp < 0 || iSp >= static_cast<int>(nSpecies) || iLv < 0) {
        SU2_MPI::Error("CONFIG ERROR: Invalid species/level index in electronic states table.", CURRENT_FUNCTION);
      }

      su2double te = 0.0, g = 0.0;
      if (c_te >= static_cast<int>(row.size()) || !TryParseDouble(row[c_te], te)) {
        SU2_MPI::Error("CONFIG ERROR: Missing characteristic electronic temperature in electronic states table.", CURRENT_FUNCTION);
      }
      if (c_g >= static_cast<int>(row.size()) || !TryParseDouble(row[c_g], g)) {
        SU2_MPI::Error("CONFIG ERROR: Missing degeneracy in electronic states table.", CURRENT_FUNCTION);
      }

      state_map[iSp][static_cast<unsigned short>(iLv)] = make_pair(te, g);
    }

    unsigned short maxEl = 0;
    for (unsigned short iSp = 0; iSp < nSpecies; ++iSp) {
      if (state_map[iSp].empty()) {
        SU2_MPI::Error("CONFIG ERROR: Every species must define at least one electronic state.", CURRENT_FUNCTION);
      }
      if (state_map[iSp].find(0) == state_map[iSp].end()) {
        SU2_MPI::Error("CONFIG ERROR: Every species must define electronic level 0.", CURRENT_FUNCTION);
      }
      nElStates[iSp] = static_cast<unsigned short>(state_map[iSp].rbegin()->first + 1);
      maxEl = max(maxEl, nElStates[iSp]);
    }

    CharElTemp.resize(nSpecies, maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies, maxEl) = su2double(0.0);

    for (unsigned short iSp = 0; iSp < nSpecies; ++iSp) {
      for (auto it = state_map[iSp].begin(); it != state_map[iSp].end(); ++it) {
        const unsigned short iLv = it->first;
        CharElTemp(iSp, iLv) = it->second.first;
        ElDegeneracy(iSp, iLv) = it->second.second;
      }
    }
  }

  /*--- Reactions table ---*/
  {
    const auto rows = read_table(reactions_file);
    const auto& header = rows[0];

    const int c_rxn = FindHeaderByTokens(header, {"reaction", "index"});
    const int c_r1 = FindHeaderByTokens(header, {"reactant", "1"});
    const int c_r2 = FindHeaderByTokens(header, {"reactant", "2"});
    const int c_r3 = FindHeaderByTokens(header, {"reactant", "3"});
    const int c_p1 = FindHeaderByTokens(header, {"product", "1"});
    const int c_p2 = FindHeaderByTokens(header, {"product", "2"});
    const int c_p3 = FindHeaderByTokens(header, {"product", "3"});
    const int c_a = FindHeaderByTokens(header, {"arrhenius_a"});
    const int c_eta = FindHeaderByTokens(header, {"arrhenius_eta"});
    const int c_theta = FindHeaderByTokens(header, {"arrhenius_theta"});
    const int c_tcf_a = FindHeaderByTokens(header, {"tcf_a"});
    const int c_tcf_b = FindHeaderByTokens(header, {"tcf_b"});
    const int c_tcb_a = FindHeaderByTokens(header, {"tcb_a"});
    const int c_tcb_b = FindHeaderByTokens(header, {"tcb_b"});
    const int c_keq_id = FindHeaderByTokens(header, {"keq", "table", "id"});

    if (c_rxn < 0 || c_r1 < 0 || c_r2 < 0 || c_r3 < 0 || c_p1 < 0 || c_p2 < 0 || c_p3 < 0 ||
        c_a < 0 || c_eta < 0 || c_theta < 0 || c_tcf_a < 0 || c_tcf_b < 0 || c_tcb_a < 0 || c_tcb_b < 0) {
      SU2_MPI::Error("CONFIG ERROR: Missing required columns in reactions table.", CURRENT_FUNCTION);
    }

    map<unsigned short, CustomReactionRow> reaction_map;

    for (size_t iRow = 1; iRow < rows.size(); ++iRow) {
      const auto& row = rows[iRow];
      if (c_rxn >= static_cast<int>(row.size())) continue;

      int rxn_idx = -1;
      if (!TryParseInt(row[c_rxn], rxn_idx)) continue;
      if (rxn_idx < 0) continue;

      CustomReactionRow reaction;
      reaction.reaction_index = static_cast<unsigned short>(rxn_idx);
      reaction.reactants[0] = ParseSpeciesSlot((c_r1 < static_cast<int>(row.size()) ? row[c_r1] : ""), nSpecies, "reactions table");
      reaction.reactants[1] = ParseSpeciesSlot((c_r2 < static_cast<int>(row.size()) ? row[c_r2] : ""), nSpecies, "reactions table");
      reaction.reactants[2] = ParseSpeciesSlot((c_r3 < static_cast<int>(row.size()) ? row[c_r3] : ""), nSpecies, "reactions table");
      reaction.products[0] = ParseSpeciesSlot((c_p1 < static_cast<int>(row.size()) ? row[c_p1] : ""), nSpecies, "reactions table");
      reaction.products[1] = ParseSpeciesSlot((c_p2 < static_cast<int>(row.size()) ? row[c_p2] : ""), nSpecies, "reactions table");
      reaction.products[2] = ParseSpeciesSlot((c_p3 < static_cast<int>(row.size()) ? row[c_p3] : ""), nSpecies, "reactions table");

      int n_react = 0, n_prod = 0;
      for (unsigned short i = 0; i < 3; ++i) {
        if (reaction.reactants[i] != nSpecies) ++n_react;
        if (reaction.products[i] != nSpecies) ++n_prod;
      }
      if (n_react == 0 || n_prod == 0) continue;

      if (c_a >= static_cast<int>(row.size()) || !TryParseDouble(row[c_a], reaction.arrhenius_a)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Arrhenius A in reactions table.", CURRENT_FUNCTION);
      }
      if (c_eta >= static_cast<int>(row.size()) || !TryParseDouble(row[c_eta], reaction.arrhenius_eta)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Arrhenius eta in reactions table.", CURRENT_FUNCTION);
      }
      if (c_theta >= static_cast<int>(row.size()) || !TryParseDouble(row[c_theta], reaction.arrhenius_theta)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Arrhenius theta in reactions table.", CURRENT_FUNCTION);
      }
      if (c_tcf_a >= static_cast<int>(row.size()) || !TryParseDouble(row[c_tcf_a], reaction.tcf_a)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Tcf_a in reactions table.", CURRENT_FUNCTION);
      }
      if (c_tcf_b >= static_cast<int>(row.size()) || !TryParseDouble(row[c_tcf_b], reaction.tcf_b)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Tcf_b in reactions table.", CURRENT_FUNCTION);
      }
      if (c_tcb_a >= static_cast<int>(row.size()) || !TryParseDouble(row[c_tcb_a], reaction.tcb_a)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Tcb_a in reactions table.", CURRENT_FUNCTION);
      }
      if (c_tcb_b >= static_cast<int>(row.size()) || !TryParseDouble(row[c_tcb_b], reaction.tcb_b)) {
        SU2_MPI::Error("CONFIG ERROR: Missing Tcb_b in reactions table.", CURRENT_FUNCTION);
      }

      if (c_keq_id >= 0 && c_keq_id < static_cast<int>(row.size())) {
        int kid = -1;
        if (TryParseInt(row[c_keq_id], kid) && kid >= 0) reaction.keq_table_id = static_cast<unsigned short>(kid);
        else reaction.keq_table_id = reaction.reaction_index;
      } else {
        reaction.keq_table_id = reaction.reaction_index;
      }

      if (reaction_map.find(reaction.reaction_index) != reaction_map.end()) {
        SU2_MPI::Error("CONFIG ERROR: Duplicate reaction_index in reactions table.", CURRENT_FUNCTION);
      }
      reaction_map[reaction.reaction_index] = reaction;
    }

    if (reaction_map.empty()) {
      SU2_MPI::Error("CONFIG ERROR: No valid reaction rows found in reactions table.", CURRENT_FUNCTION);
    }

    nReactions = static_cast<unsigned short>(reaction_map.rbegin()->first + 1);
    if (reaction_map.size() != nReactions) {
      SU2_MPI::Error("CONFIG ERROR: Reactions must use contiguous reaction_index values starting from 0.", CURRENT_FUNCTION);
    }

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.assign(nReactions, 0.0);
    ArrheniusEta.assign(nReactions, 0.0);
    ArrheniusTheta.assign(nReactions, 0.0);
    Tcf_a.assign(nReactions, 0.0);
    Tcf_b.assign(nReactions, 0.0);
    Tcb_a.assign(nReactions, 0.0);
    Tcb_b.assign(nReactions, 0.0);
    KeqTableIndex.assign(nReactions, 0);

    for (unsigned short iRxn = 0; iRxn < nReactions; ++iRxn) {
      for (unsigned short i = 0; i < 3; ++i) {
        Reactions(iRxn,0,i) = nSpecies;
        Reactions(iRxn,1,i) = nSpecies;
      }
      const auto& reaction = reaction_map[iRxn];
      for (unsigned short i = 0; i < 3; ++i) {
        Reactions(iRxn,0,i) = reaction.reactants[i];
        Reactions(iRxn,1,i) = reaction.products[i];
      }

      ArrheniusCoefficient[iRxn] = reaction.arrhenius_a;
      ArrheniusEta[iRxn] = reaction.arrhenius_eta;
      ArrheniusTheta[iRxn] = reaction.arrhenius_theta;
      Tcf_a[iRxn] = reaction.tcf_a;
      Tcf_b[iRxn] = reaction.tcf_b;
      Tcb_a[iRxn] = reaction.tcb_a;
      Tcb_b[iRxn] = reaction.tcb_b;
      KeqTableIndex[iRxn] = reaction.keq_table_id;
    }
  }

  /*--- Keq coefficient table: allow a single coefficient set per reaction (no density segmentation) ---*/
  {
    const auto rows = read_table(keq_file);
    const auto& header = rows[0];

    const int c_keq_id = FindHeaderByTokens(header, {"keq", "table", "id"});
    const int c_rxn = FindHeaderByTokens(header, {"reaction", "index"});
    const int c_a0 = FindHeaderByTokens(header, {"a0"});
    const int c_a1 = FindHeaderByTokens(header, {"a1"});
    const int c_a2 = FindHeaderByTokens(header, {"a2"});
    const int c_a3 = FindHeaderByTokens(header, {"a3"});
    const int c_a4 = FindHeaderByTokens(header, {"a4"});

    std::unordered_map<unsigned short, std::array<su2double,5>> coeff_by_key;

    auto set_coeff = [&coeff_by_key](unsigned short key, const std::array<su2double,5>& coeff) {
      auto found = coeff_by_key.find(key);
      if (found == coeff_by_key.end()) {
        coeff_by_key[key] = coeff;
        return;
      }
      for (unsigned short i = 0; i < 5; ++i) {
        if (fabs(found->second[i] - coeff[i]) > 1e-10) {
          SU2_MPI::Error("CONFIG ERROR: Multiple different Keq coefficient sets provided for same reaction/table id."
                         " For no-density-segmentation mode, only one set is allowed.", CURRENT_FUNCTION);
        }
      }
    };

    for (size_t iRow = 1; iRow < rows.size(); ++iRow) {
      const auto& row = rows[iRow];

      bool parsed = false;

      if (c_a0 >= 0 && c_a1 >= 0 && c_a2 >= 0 && c_a3 >= 0 && c_a4 >= 0 &&
          c_a4 < static_cast<int>(row.size())) {
        su2double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
        if (TryParseDouble(row[c_a0], a0) && TryParseDouble(row[c_a1], a1) &&
            TryParseDouble(row[c_a2], a2) && TryParseDouble(row[c_a3], a3) &&
            TryParseDouble(row[c_a4], a4)) {
          int key = -1;
          if (c_keq_id >= 0 && c_keq_id < static_cast<int>(row.size()) && TryParseInt(row[c_keq_id], key) && key >= 0) {
            parsed = true;
          } else if (c_rxn >= 0 && c_rxn < static_cast<int>(row.size()) && TryParseInt(row[c_rxn], key) && key >= 0) {
            parsed = true;
          }

          if (parsed) {
            set_coeff(static_cast<unsigned short>(key), std::array<su2double,5>{{a0,a1,a2,a3,a4}});
          }
        }
      }

      if (!parsed && row.size() >= 7) {
        int key = -1;
        su2double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
        if (TryParseInt(row[0], key) && key >= 0 &&
            TryParseDouble(row[2], a0) && TryParseDouble(row[3], a1) &&
            TryParseDouble(row[4], a2) && TryParseDouble(row[5], a3) &&
            TryParseDouble(row[6], a4)) {
          set_coeff(static_cast<unsigned short>(key), std::array<su2double,5>{{a0,a1,a2,a3,a4}});
        }
      }
    }

    if (coeff_by_key.empty()) {
      SU2_MPI::Error("CONFIG ERROR: No valid Keq coefficients found in Keq table.", CURRENT_FUNCTION);
    }

    CustomKeqCoeff.assign(nReactions, std::array<su2double,5>{{0.0,0.0,0.0,0.0,0.0}});
    for (unsigned short iRxn = 0; iRxn < nReactions; ++iRxn) {
      const unsigned short table_key = KeqTableIndex[iRxn];
      auto found = coeff_by_key.find(table_key);
      if (found == coeff_by_key.end()) {
        found = coeff_by_key.find(iRxn);
      }
      if (found == coeff_by_key.end()) {
        SU2_MPI::Error("CONFIG ERROR: Missing Keq coefficient set for reaction index " + std::to_string(iRxn), CURRENT_FUNCTION);
      }
      CustomKeqCoeff[iRxn] = found->second;
    }
  }

  /*--- Transport table ---*/
  if (viscous) {
    const auto rows = read_table(transport_file);
    const auto& header = rows[0];

    const int c_model = FindHeaderByTokens(header, {"model_type"});
    const int c_si = FindHeaderByTokens(header, {"species_i"});
    const int c_sj = FindHeaderByTokens(header, {"species_j"});
    const int c_o110 = FindHeaderByTokens(header, {"omega11", "0"});
    const int c_o111 = FindHeaderByTokens(header, {"omega11", "1"});
    const int c_o112 = FindHeaderByTokens(header, {"omega11", "2"});
    const int c_o113 = FindHeaderByTokens(header, {"omega11", "3"});
    const int c_o220 = FindHeaderByTokens(header, {"omega22", "0"});
    const int c_o221 = FindHeaderByTokens(header, {"omega22", "1"});
    const int c_o222 = FindHeaderByTokens(header, {"omega22", "2"});
    const int c_o223 = FindHeaderByTokens(header, {"omega22", "3"});
    const int c_ba = FindHeaderByTokens(header, {"blottner_a"});
    const int c_bb = FindHeaderByTokens(header, {"blottner_b"});
    const int c_bc = FindHeaderByTokens(header, {"blottner_c"});
    const int c_mu = FindHeaderByTokens(header, {"mu_ref"});
    const int c_k = FindHeaderByTokens(header, {"k_ref"});
    const int c_sm = FindHeaderByTokens(header, {"sm_ref"});
    const int c_sk = FindHeaderByTokens(header, {"sk_ref"});

    if (c_model < 0) {
      SU2_MPI::Error("CONFIG ERROR: Missing model_type column in transport table.", CURRENT_FUNCTION);
    }

    vector<bool> blottner_set(nSpecies, false);
    vector<vector<bool>> omega11_set(nSpecies, vector<bool>(nSpecies, false));
    bool has_sutherland = false;

    for (size_t iRow = 1; iRow < rows.size(); ++iRow) {
      const auto& row = rows[iRow];
      if (c_model >= static_cast<int>(row.size())) continue;
      const string model = CanonicalHeader(row[c_model]);

      if (model.find("sutherland") != string::npos) {
        su2double mu_ref_val = 0.0, k_ref_val = 0.0, sm_ref_val = 0.0, sk_ref_val = 0.0;
        if (c_mu < static_cast<int>(row.size()) && c_k < static_cast<int>(row.size()) &&
            c_sm < static_cast<int>(row.size()) && c_sk < static_cast<int>(row.size()) &&
            TryParseDouble(row[c_mu], mu_ref_val) && TryParseDouble(row[c_k], k_ref_val) &&
            TryParseDouble(row[c_sm], sm_ref_val) && TryParseDouble(row[c_sk], sk_ref_val)) {
          mu_ref[0] = mu_ref_val;
          k_ref[0] = k_ref_val;
          Sm_ref[0] = sm_ref_val;
          Sk_ref[0] = sk_ref_val;
          has_sutherland = true;
        }
        continue;
      }

      int si = -1, sj = -1;
      if (c_si < 0 || c_sj < 0 || c_si >= static_cast<int>(row.size()) || c_sj >= static_cast<int>(row.size())) continue;
      if (!TryParseInt(row[c_si], si) || !TryParseInt(row[c_sj], sj)) continue;
      if (si < 0 || sj < 0 || si >= static_cast<int>(nSpecies) || sj >= static_cast<int>(nSpecies)) {
        SU2_MPI::Error("CONFIG ERROR: species_i/species_j out of range in transport table.", CURRENT_FUNCTION);
      }

      su2double o110 = 0.0, o111 = 0.0, o112 = 0.0, o113 = 0.0;
      su2double o220 = 0.0, o221 = 0.0, o222 = 0.0, o223 = 0.0;
      if (c_o110 >= 0 && c_o111 >= 0 && c_o112 >= 0 && c_o113 >= 0 &&
          c_o113 < static_cast<int>(row.size()) &&
          TryParseDouble(row[c_o110], o110) && TryParseDouble(row[c_o111], o111) &&
          TryParseDouble(row[c_o112], o112) && TryParseDouble(row[c_o113], o113)) {
        Omega11(si,sj,0) = o110; Omega11(si,sj,1) = o111; Omega11(si,sj,2) = o112; Omega11(si,sj,3) = o113;
        Omega11(sj,si,0) = o110; Omega11(sj,si,1) = o111; Omega11(sj,si,2) = o112; Omega11(sj,si,3) = o113;
        omega11_set[si][sj] = true;
        omega11_set[sj][si] = true;
      }

      if (c_o220 >= 0 && c_o221 >= 0 && c_o222 >= 0 && c_o223 >= 0 &&
          c_o223 < static_cast<int>(row.size()) &&
          TryParseDouble(row[c_o220], o220) && TryParseDouble(row[c_o221], o221) &&
          TryParseDouble(row[c_o222], o222) && TryParseDouble(row[c_o223], o223)) {
        Omega22(si,sj,0) = o220; Omega22(si,sj,1) = o221; Omega22(si,sj,2) = o222; Omega22(si,sj,3) = o223;
        Omega22(sj,si,0) = o220; Omega22(sj,si,1) = o221; Omega22(sj,si,2) = o222; Omega22(sj,si,3) = o223;
      }

      su2double ba = 0.0, bb = 0.0, bc = 0.0;
      if (c_ba >= 0 && c_bb >= 0 && c_bc >= 0 && c_bc < static_cast<int>(row.size()) &&
          TryParseDouble(row[c_ba], ba) && TryParseDouble(row[c_bb], bb) && TryParseDouble(row[c_bc], bc)) {
        Blottner(si,0) = ba;
        Blottner(si,1) = bb;
        Blottner(si,2) = bc;
        blottner_set[si] = true;
      }
    }

    if (Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND && !has_sutherland) {
      SU2_MPI::Error("CONFIG ERROR: SUTHERLAND transport selected but no valid Sutherland coefficients found in transport table.", CURRENT_FUNCTION);
    }

    if (Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE) {
      for (unsigned short iSp = 0; iSp < nSpecies; ++iSp) {
        if (!blottner_set[iSp]) {
          SU2_MPI::Error("CONFIG ERROR: Missing Blottner coefficients for one or more species in transport table.", CURRENT_FUNCTION);
        }
      }

      for (unsigned short iSp = 0; iSp < nSpecies; ++iSp) {
        for (unsigned short jSp = iSp; jSp < nSpecies; ++jSp) {
          if (!omega11_set[iSp][jSp]) {
            SU2_MPI::Error("CONFIG ERROR: Missing Omega11 coefficients for species pair (" +
                               std::to_string(iSp) + "," + std::to_string(jSp) + ") in transport table.",
                           CURRENT_FUNCTION);
          }
        }
      }
    }
  }
}

CSU2TCLib::CSU2TCLib(const CConfig* config, unsigned short val_nDim, bool viscous): CNEMOGas(config, val_nDim){

  unsigned short maxEl = 0;
  su2double mf = 0.0;

  const auto MassFrac_Freestream = config->GetGas_Composition();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mf += MassFrac_Freestream[iSpecies];

  /*--- Allocate vectors for gas properties ---*/
  nElStates.resize(nSpecies,0);
  CharVibTemp.resize(nSpecies,0.0);
  RotationModes.resize(nSpecies,0.0);
  Diss.resize(nSpecies,0.0);
  A.resize(5,0.0);
  Omega11.resize(nSpecies,nSpecies,4,0.0);
  Omega22.resize(nSpecies,nSpecies,4,0.0);
  RxnConstantTable.resize(6,5) = su2double(0.0);
  CatRecombTable.resize(nSpecies,2) = 0;
  Blottner.resize(nSpecies,3)  = su2double(0.0);
  taus.resize(nSpecies,0.0);
  eve_eq.resize(nSpecies,0.0);
  eve.resize(nSpecies,0.0);

  if (viscous) {
    MolarFracWBE.resize(nSpecies,0.0);
    phis.resize(nSpecies,0.0);
    mus.resize(nSpecies,0.0);
  }

  custom_chemistry_tables = config->GetSU2_NONEQ_CustomChemistry();
  if (custom_chemistry_tables) {
    LoadCustomChemistryTables(config, viscous);

    ionization = false;
    nHeavy = nSpecies;
    nEl = 0;

    if (nReactions == 0) {
      SU2_MPI::Error("CONFIG ERROR: Custom chemistry mode requires at least one reaction.", CURRENT_FUNCTION);
    }

    return;
  }

  if (gas_model =="ARGON"){
    if (nSpecies != 1) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.667;
    nReactions  = 0;

    // Molar mass [kg/kmol]
    MolarMass[0] = 39.948;
    // Rotational modes of energy storage
    RotationModes[0] = 0.0;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;

    Enthalpy_Formation[0] = 0.0;
    Ref_Temperature[0] = 0.0;
    nElStates[0] = 7;

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- AR: Blottner coefficients. ---*/
    Blottner(0,0) = 3.83444322E-03;   Blottner(0,1) = 6.74718764E-01;   Blottner(0,2) = -1.24290388E+01;

    /*--- AR: 7 states ---*/
    CharElTemp(0,0) = 0.000000000000000E+00;
    CharElTemp(0,1) = 1.611135736988230E+05;
    CharElTemp(0,2) = 1.625833076870950E+05;
    CharElTemp(0,3) = 1.636126382960720E+05;
    CharElTemp(0,4) = 1.642329518358000E+05;
    CharElTemp(0,5) = 1.649426852542080E+05;
    CharElTemp(0,6) = 1.653517702884570E+05;
    ElDegeneracy(0,0) = 1;
    ElDegeneracy(0,1) = 9;
    ElDegeneracy(0,2) = 21;
    ElDegeneracy(0,3) = 7;
    ElDegeneracy(0,4) = 3;
    ElDegeneracy(0,5) = 5;
    ElDegeneracy(0,6) = 15;

    /*--- Catalytic wall table ---*/
    // Creation/Destruction (+1/-1), Index of monoatomic reactants.
    // Argon not used.
    CatRecombTable(0,0) = 0; CatRecombTable(0,1) = 0;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      mu_ref[0] = 2.125E-5;
      k_ref[0] = 0.0163;
      Sm_ref[0] = 114.0;
      Sk_ref[0] = 170;
    }

  } else if (gas_model == "N2"){
    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 2) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 2;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 14.0067;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 0.0;
    // Formation enthalpy: (JANAF values [KJ/Kmol])
    // J/kg - from Scalabrin
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 3.36E7;   //N
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    // Blottner viscosity coefficients
    // A                       // B                       // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) = 3.18E-1;   Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 1.16E-2;   Blottner(1,1) = 6.03E-1;   Blottner(1,2) = -1.24E1;  // N
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 3;                     // N
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- Assign values to data structures ---*/
    // N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // N: 3 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 2.766469645581980E+04;
    CharElTemp(1,2) = 4.149309313560210E+04;
    ElDegeneracy(1,0) = 4;
    ElDegeneracy(1,1) = 10;
    ElDegeneracy(1,2) = 6;
    /*--- Set Arrhenius coefficients for chemical reactions ---*/
    // Note: Data lists coefficients in (cm^3/mol-s) units, need to convert
    //       to (m^3/kmol-s) to be consistent with the rest of the code
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 3.0E22;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    // Characteristic temperature
    ArrheniusTheta[0] = 113200.0;
    ArrheniusTheta[1] = 113200.0;
    /*--- Set reaction maps ---*/
    // N2 + N2 -> 2N + N2
    Reactions(0,0,0)=0;   Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;
    Reactions(0,1,0)=1;   Reactions(0,1,1)=1;   Reactions(0,1,2) =0;
    // N2 + N -> 2N + N
    Reactions(1,0,0)=0;   Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;
    Reactions(1,1,0)=1;   Reactions(1,1,1)=1;   Reactions(1,1,2)=1;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0] = 0.5; Tcf_b[0] = 0.5; Tcb_a[0] = 1;  Tcb_b[0] = 0;
    Tcf_a[1] = 0.5; Tcf_b[1] = 0.5; Tcb_a[1] = 1;  Tcb_b[1] = 0;

    /*--- Dissociation potential [KJ/kg] ---*/
    Diss[0] = 3.36E4;
    Diss[1] = 0.0;

    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3
    Omega11(0,0,0) = -6.0614558E-03;  Omega11(0,0,1) = 1.2689102E-01;   Omega11(0,0,2) = -1.0616948E+00;  Omega11(0,0,3) = 8.0955466E+02;
    Omega11(0,1,0) = -1.0796249E-02;  Omega11(0,1,1) = 2.2656509E-01;   Omega11(0,1,2) = -1.7910602E+00;  Omega11(0,1,3) = 4.0455218E+03;
    Omega11(1,0,0) = -1.0796249E-02;  Omega11(1,0,1) = 2.2656509E-01;   Omega11(1,0,2) = -1.7910602E+00;  Omega11(1,0,3) = 4.0455218E+03;
    Omega11(1,1,0) = -9.6083779E-03;  Omega11(1,1,1) = 2.0938971E-01;   Omega11(1,1,2) = -1.7386904E+00;  Omega11(1,1,3) = 3.3587983E+03;
    Omega22(0,0,0) = -7.6303990E-03;  Omega22(0,0,1) = 1.6878089E-01;   Omega22(0,0,2) = -1.4004234E+00;  Omega22(0,0,3) = 2.1427708E+03;
    Omega22(0,1,0) = -8.3493693E-03;  Omega22(0,1,1) = 1.7808911E-01;   Omega22(0,1,2) = -1.4466155E+00;  Omega22(0,1,3) = 1.9324210E+03;
    Omega22(1,0,0) = -8.3493693E-03;  Omega22(1,0,1) = 1.7808911E-01;   Omega22(1,0,2) = -1.4466155E+00;  Omega22(1,0,3) = 1.9324210E+03;
    Omega22(1,1,0) = -7.7439615E-03;  Omega22(1,1,1) = 1.7129007E-01;   Omega22(1,1,2) = -1.4809088E+00;  Omega22(1,1,3) = 2.1284951E+03;

    /*--- Catalytic wall table ---*/
    // Creation/Destruction (+1/-1), Index of monoatomic reactants.
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  1; CatRecombTable(0,1) = 1;
    CatRecombTable(1,0) = -1; CatRecombTable(1,1) = 1;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0242;
      mu_ref[0] = 1.663E-5;
      Sm_ref[0] = 107.0;
      Sk_ref[0] = 150.0;
    }

  } else if (gas_model == "AIR-5"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 5) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition",CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 17;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 2.0;
    RotationModes[2] = 2.0;
    RotationModes[3] = 0.0;
    RotationModes[4] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 2.0*15.9994;
    MolarMass[2] = 14.0067+15.9994;
    MolarMass[3] = 14.0067;
    MolarMass[4] = 15.9994;
    //Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 2239.0;
    CharVibTemp[2] = 2817.0;
    CharVibTemp[3] = 0.0;
    CharVibTemp[4] = 0.0;
    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 0.0;      //O2
    Enthalpy_Formation[2] = 3.0E6;    //NO
    Enthalpy_Formation[3] = 3.36E7;   //N
    Enthalpy_Formation[4] = 1.54E7;   //O
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) =  3.18E-1;  Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 4.49E-2;   Blottner(1,1) = -8.26E-2;  Blottner(1,2) = -9.20E0;  // O2
    Blottner(2,0) = 4.36E-2;   Blottner(2,1) = -3.36E-2;  Blottner(2,2) = -9.58E0;  // NO
    Blottner(3,0) = 1.16E-2;   Blottner(3,1) =  6.03E-1;  Blottner(3,2) = -1.24E1;  // N
    Blottner(4,0) = 2.03E-2;   Blottner(4,1) =  4.29E-1;  Blottner(4,2) = -1.16E1;  // O
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 7;                     // O2
    nElStates[2] = 16;                    // NO
    nElStates[3] = 3;                     // N
    nElStates[4] = 5;                     // O
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);
    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    //N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // O2: 7 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 1.139156019700800E+04;
    CharElTemp(1,2) = 1.898473947826400E+04;
    CharElTemp(1,3) = 4.755973576639200E+04;
    CharElTemp(1,4) = 4.991242097343200E+04;
    CharElTemp(1,5) = 5.092268575561600E+04;
    CharElTemp(1,6) = 7.189863255967200E+04;
    ElDegeneracy(1,0) = 3;
    ElDegeneracy(1,1) = 2;
    ElDegeneracy(1,2) = 1;
    ElDegeneracy(1,3) = 1;
    ElDegeneracy(1,4) = 6;
    ElDegeneracy(1,5) = 3;
    ElDegeneracy(1,6) = 3;
    // NO: 16 states
    CharElTemp(2,0)  = 0.000000000000000E+00;
    CharElTemp(2,1)  = 5.467345760000000E+04;
    CharElTemp(2,2)  = 6.317139627802400E+04;
    CharElTemp(2,3)  = 6.599450342445600E+04;
    CharElTemp(2,4)  = 6.906120960000000E+04;
    CharElTemp(2,5)  = 7.049998480000000E+04;
    CharElTemp(2,6)  = 7.491055017560000E+04;
    CharElTemp(2,7)  = 7.628875293968000E+04;
    CharElTemp(2,8)  = 8.676188537552000E+04;
    CharElTemp(2,9)  = 8.714431182368000E+04;
    CharElTemp(2,10) = 8.886077063728000E+04;
    CharElTemp(2,11) = 8.981755614528000E+04;
    CharElTemp(2,12) = 8.988445919208000E+04;
    CharElTemp(2,13) = 9.042702132000000E+04;
    CharElTemp(2,14) = 9.064283760000000E+04;
    CharElTemp(2,15) = 9.111763341600000E+04;
    ElDegeneracy(2,0)  = 4;
    ElDegeneracy(2,1)  = 8;
    ElDegeneracy(2,2)  = 2;
    ElDegeneracy(2,3)  = 4;
    ElDegeneracy(2,4)  = 4;
    ElDegeneracy(2,5)  = 4;
    ElDegeneracy(2,6)  = 4;
    ElDegeneracy(2,7)  = 2;
    ElDegeneracy(2,8)  = 4;
    ElDegeneracy(2,9)  = 2;
    ElDegeneracy(2,10) = 4;
    ElDegeneracy(2,11) = 4;
    ElDegeneracy(2,12) = 2;
    ElDegeneracy(2,13) = 2;
    ElDegeneracy(2,14) = 2;
    ElDegeneracy(2,15) = 4;
    // N: 3 states
    CharElTemp(3,0) = 0.000000000000000E+00;
    CharElTemp(3,1) = 2.766469645581980E+04;
    CharElTemp(3,2) = 4.149309313560210E+04;
    ElDegeneracy(3,0)= 4;
    ElDegeneracy(3,1)= 10;
    ElDegeneracy(3,2)= 6;
    // O: 5 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.277077570280000E+02;
    CharElTemp(4,2) = 3.265688785704000E+02;
    CharElTemp(4,3) = 2.283028632262240E+04;
    CharElTemp(4,4) = 4.861993036434160E+04;
    ElDegeneracy(4,0) = 5;
    ElDegeneracy(4,1) = 3;
    ElDegeneracy(4,2) = 1;
    ElDegeneracy(4,3) = 5;
    ElDegeneracy(4,4) = 1;
    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions(0,0,0)=0;    Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;    Reactions(0,1,0)=3;   Reactions(0,1,1)=3;   Reactions(0,1,2) =0;
    Reactions(1,0,0)=0;    Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;    Reactions(1,1,0)=3;   Reactions(1,1,1)=3;   Reactions(1,1,2) =1;
    Reactions(2,0,0)=0;    Reactions(2,0,1)=2;   Reactions(2,0,2)=nSpecies;    Reactions(2,1,0)=3;   Reactions(2,1,1)=3;   Reactions(2,1,2) =2;
    Reactions(3,0,0)=0;    Reactions(3,0,1)=3;   Reactions(3,0,2)=nSpecies;    Reactions(3,1,0)=3;   Reactions(3,1,1)=3;   Reactions(3,1,2) =3;
    Reactions(4,0,0)=0;    Reactions(4,0,1)=4;   Reactions(4,0,2)=nSpecies;    Reactions(4,1,0)=3;   Reactions(4,1,1)=3;   Reactions(4,1,2) =4;
    // O2 dissociation
    Reactions(5,0,0)=1;    Reactions(5,0,1)=0;   Reactions(5,0,2)=nSpecies;    Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2) =0;
    Reactions(6,0,0)=1;    Reactions(6,0,1)=1;   Reactions(6,0,2)=nSpecies;    Reactions(6,1,0)=4;   Reactions(6,1,1)=4;   Reactions(6,1,2) =1;
    Reactions(7,0,0)=1;    Reactions(7,0,1)=2;   Reactions(7,0,2)=nSpecies;    Reactions(7,1,0)=4;   Reactions(7,1,1)=4;   Reactions(7,1,2) =2;
    Reactions(8,0,0)=1;    Reactions(8,0,1)=3;   Reactions(8,0,2)=nSpecies;    Reactions(8,1,0)=4;   Reactions(8,1,1)=4;   Reactions(8,1,2) =3;
    Reactions(9,0,0)=1;    Reactions(9,0,1)=4;   Reactions(9,0,2)=nSpecies;    Reactions(9,1,0)=4;   Reactions(9,1,1)=4;   Reactions(9,1,2) =4;
    // NO dissociation
    Reactions(10,0,0)=2;   Reactions(10,0,1)=0;  Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=3;  Reactions(10,1,1)=4;    Reactions(10,1,2) =0;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=1;  Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=3;  Reactions(11,1,1)=4;    Reactions(11,1,2) =1;
    Reactions(12,0,0)=2;   Reactions(12,0,1)=2;  Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=3;  Reactions(12,1,1)=4;    Reactions(12,1,2) =2;
    Reactions(13,0,0)=2;   Reactions(13,0,1)=3;  Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=3;  Reactions(13,1,1)=4;    Reactions(13,1,2) =3;
    Reactions(14,0,0)=2;   Reactions(14,0,1)=4;  Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=3;  Reactions(14,1,1)=4;    Reactions(14,1,2) =4;
    // N2 + O -> NO + N
    Reactions(15,0,0)=0;   Reactions(15,0,1)=4;  Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=2;  Reactions(15,1,1)=3;    Reactions(15,1,2)= nSpecies;
    // NO + O -> O2 + N
    Reactions(16,0,0)=2;   Reactions(16,0,1)=4;  Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=1;  Reactions(16,1,1)=3;    Reactions(16,1,2)= nSpecies;
    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 2.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 1.0E22;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 5.0E15;
    ArrheniusCoefficient[11] = 5.0E15;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 1.1E17;
    ArrheniusCoefficient[14] = 1.1E17;
    ArrheniusCoefficient[15] = 6.4E17;
    ArrheniusCoefficient[16] = 8.4E12;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.50;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = 0.0;
    ArrheniusEta[11] = 0.0;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = -1.0;
    ArrheniusEta[16] = 0.0;
    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 59500.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10] = 75500.0;
    ArrheniusTheta[11] = 75500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 38400.0;
    ArrheniusTheta[16] = 19450.0;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;  Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;  Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;  Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;  Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;  Tcb_b[4] = 0;
    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;  Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;  Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;  Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;  Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;  Tcb_b[9] = 0;
    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;  Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;  Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;  Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;  Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;  Tcb_b[14] = 0;
    Tcf_a[15] = 1.0; Tcf_b[15] = 0.0; Tcb_a[15] = 1;  Tcb_b[15] = 0;
    Tcf_a[16] = 1.0; Tcf_b[16] = 0.0; Tcb_a[16] = 1;  Tcb_b[16] = 0;
    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3
    // Omega^(1,1) ----------------------
    //N2
    Omega11(0,0,0) = -6.0614558E-03;  Omega11(0,0,1) = 1.2689102E-01;   Omega11(0,0,2) = -1.0616948E+00;  Omega11(0,0,3) = 8.0955466E+02;
    Omega11(0,1,0) = -3.7959091E-03;  Omega11(0,1,1) = 9.5708295E-02;   Omega11(0,1,2) = -1.0070611E+00;  Omega11(0,1,3) = 8.9392313E+02;
    Omega11(0,2,0) = -1.9295666E-03;  Omega11(0,2,1) = 2.7995735E-02;   Omega11(0,2,2) = -3.1588514E-01;  Omega11(0,2,3) = 1.2880734E+02;
    Omega11(0,3,0) = -1.0796249E-02;  Omega11(0,3,1) = 2.2656509E-01;   Omega11(0,3,2) = -1.7910602E+00;  Omega11(0,3,3) = 4.0455218E+03;
    Omega11(0,4,0) = -2.7244269E-03;  Omega11(0,4,1) = 6.9587171E-02;   Omega11(0,4,2) = -7.9538667E-01;  Omega11(0,4,3) = 4.0673730E+02;
    //O2
    Omega11(1,0,0) = -3.7959091E-03;  Omega11(1,0,1) = 9.5708295E-02;   Omega11(1,0,2) = -1.0070611E+00;  Omega11(1,0,3) = 8.9392313E+02;
    Omega11(1,1,0) = -8.0682650E-04;  Omega11(1,1,1) = 1.6602480E-02;   Omega11(1,1,2) = -3.1472774E-01;  Omega11(1,1,3) = 1.4116458E+02;
    Omega11(1,2,0) = -6.4433840E-04;  Omega11(1,2,1) = 8.5378580E-03;   Omega11(1,2,2) = -2.3225102E-01;  Omega11(1,2,3) = 1.1371608E+02;
    Omega11(1,3,0) = -1.1453028E-03;  Omega11(1,3,1) = 1.2654140E-02;   Omega11(1,3,2) = -2.2435218E-01;  Omega11(1,3,3) = 7.7201588E+01;
    Omega11(1,4,0) = -4.8405803E-03;  Omega11(1,4,1) = 1.0297688E-01;   Omega11(1,4,2) = -9.6876576E-01;  Omega11(1,4,3) = 6.1629812E+02;
    //NO
    Omega11(2,0,0) = -1.9295666E-03;  Omega11(2,0,1) = 2.7995735E-02;   Omega11(2,0,2) = -3.1588514E-01;  Omega11(2,0,3) = 1.2880734E+02;
    Omega11(2,1,0) = -6.4433840E-04;  Omega11(2,1,1) = 8.5378580E-03;   Omega11(2,1,2) = -2.3225102E-01;  Omega11(2,1,3) = 1.1371608E+02;
    Omega11(2,2,0) = -0.0000000E+00;  Omega11(2,2,1) = -1.1056066E-02;  Omega11(2,2,2) = -5.9216250E-02;  Omega11(2,2,3) = 7.2542367E+01;
    Omega11(2,3,0) = -1.5770918E-03;  Omega11(2,3,1) = 1.9578381E-02;   Omega11(2,3,2) = -2.7873624E-01;  Omega11(2,3,3) = 9.9547944E+01;
    Omega11(2,4,0) = -1.0885815E-03;  Omega11(2,4,1) = 1.1883688E-02;   Omega11(2,4,2) = -2.1844909E-01;  Omega11(2,4,3) = 7.5512560E+01;
    //N
    Omega11(3,0,0) = -1.0796249E-02;  Omega11(3,0,1) = 2.2656509E-01;   Omega11(3,0,2) = -1.7910602E+00;  Omega11(3,0,3) = 4.0455218E+03;
    Omega11(3,1,0) = -1.1453028E-03;  Omega11(3,1,1) = 1.2654140E-02;   Omega11(3,1,2) = -2.2435218E-01;  Omega11(3,1,3) = 7.7201588E+01;
    Omega11(3,2,0) = -1.5770918E-03;  Omega11(3,2,1) = 1.9578381E-02;   Omega11(3,2,2) = -2.7873624E-01;  Omega11(3,2,3) = 9.9547944E+01;
    Omega11(3,3,0) = -9.6083779E-03;  Omega11(3,3,1) = 2.0938971E-01;   Omega11(3,3,2) = -1.7386904E+00;  Omega11(3,3,3) = 3.3587983E+03;
    Omega11(3,4,0) = -7.8147689E-03;  Omega11(3,4,1) = 1.6792705E-01;   Omega11(3,4,2) = -1.4308628E+00;  Omega11(3,4,3) = 1.6628859E+03;
    //O
    Omega11(4,0,0) = -2.7244269E-03;  Omega11(4,0,1) = 6.9587171E-02;   Omega11(4,0,2) = -7.9538667E-01;  Omega11(4,0,3) = 4.0673730E+02;
    Omega11(4,1,0) = -4.8405803E-03;  Omega11(4,1,1) = 1.0297688E-01;   Omega11(4,1,2) = -9.6876576E-01;  Omega11(4,1,3) = 6.1629812E+02;
    Omega11(4,2,0) = -1.0885815E-03;  Omega11(4,2,1) = 1.1883688E-02;   Omega11(4,2,2) = -2.1844909E-01;  Omega11(4,2,3) = 7.5512560E+01;
    Omega11(4,3,0) = -7.8147689E-03;  Omega11(4,3,1) = 1.6792705E-01;   Omega11(4,3,2) = -1.4308628E+00;  Omega11(4,3,3) = 1.6628859E+03;
    Omega11(4,4,0) = -6.4040535E-03;  Omega11(4,4,1) = 1.4629949E-01;   Omega11(4,4,2) = -1.3892121E+00;  Omega11(4,4,3) = 2.0903441E+03;
 
    // Omega^(2,2) ----------------------
    //N2
    Omega22(0,0,0) = -7.6303990E-03;  Omega22(0,0,1) = 1.6878089E-01;   Omega22(0,0,2) = -1.4004234E+00;  Omega22(0,0,3) = 2.1427708E+03;
    Omega22(0,1,0) = -8.0457321E-03;  Omega22(0,1,1) = 1.9228905E-01;   Omega22(0,1,2) = -1.7102854E+00;  Omega22(0,1,3) = 5.2213857E+03;
    Omega22(0,2,0) = -6.8237776E-03;  Omega22(0,2,1) = 1.4360616E-01;   Omega22(0,2,2) = -1.1922240E+00;  Omega22(0,2,3) = 1.2433086E+03;
    Omega22(0,3,0) = -8.3493693E-03;  Omega22(0,3,1) = 1.7808911E-01;   Omega22(0,3,2) = -1.4466155E+00;  Omega22(0,3,3) = 1.9324210E+03;
    Omega22(0,4,0) = -8.3110691E-03;  Omega22(0,4,1) = 1.9617877E-01;   Omega22(0,4,2) = -1.7205427E+00;  Omega22(0,4,3) = 4.0812829E+03;
    //O2
    Omega22(1,0,0) = -8.0457321E-03;  Omega22(1,0,1) = 1.9228905E-01;   Omega22(1,0,2) = -1.7102854E+00;  Omega22(1,0,3) = 5.2213857E+03;
    Omega22(1,1,0) = -6.2931612E-03;  Omega22(1,1,1) = 1.4624645E-01;   Omega22(1,1,2) = -1.3006927E+00;  Omega22(1,1,3) = 1.8066892E+03;
    Omega22(1,2,0) = -6.8508672E-03;  Omega22(1,2,1) = 1.5524564E-01;   Omega22(1,2,2) = -1.3479583E+00;  Omega22(1,2,3) = 2.0037890E+03;
    Omega22(1,3,0) = -1.0608832E-03;  Omega22(1,3,1) = 1.1782595E-02;   Omega22(1,3,2) = -2.1246301E-01;  Omega22(1,3,3) = 8.4561598E+01;
    Omega22(1,4,0) = -3.7969686E-03;  Omega22(1,4,1) = 7.6789981E-02;   Omega22(1,4,2) = -7.3056809E-01;  Omega22(1,4,3) = 3.3958171E+02;
    //NO
    Omega22(2,0,0) = -6.8237776E-03;  Omega22(2,0,1) = 1.4360616E-01;   Omega22(2,0,2) = -1.1922240E+00;  Omega22(2,0,3) = 1.2433086E+03;
    Omega22(2,1,0) = -6.8508672E-03;  Omega22(2,1,1) = 1.5524564E-01;   Omega22(2,1,2) = -1.3479583E+00;  Omega22(2,1,3) = 2.0037890E+03;
    Omega22(2,2,0) = -7.4942466E-03;  Omega22(2,2,1) = 1.6626193E-01;   Omega22(2,2,2) = -1.4107027E+00;  Omega22(2,2,3) = 2.3097604E+03;
    Omega22(2,3,0) = -1.4719259E-03;  Omega22(2,3,1) = 1.8446968E-02;   Omega22(2,3,2) = -2.6460411E-01;  Omega22(2,3,3) = 1.0911124E+02;
    Omega22(2,4,0) = -1.0066279E-03;  Omega22(2,4,1) = 1.1029264E-02;   Omega22(2,4,2) = -2.0671266E-01;  Omega22(2,4,3) = 8.2644384E+01;
    //N
    Omega22(3,0,0) = -8.3493693E-03;  Omega22(3,0,1) = 1.7808911E-01;   Omega22(3,0,2) = -1.4466155E+00;  Omega22(3,0,3) = 1.9324210E+03;
    Omega22(3,1,0) = -1.0608832E-03;  Omega22(3,1,1) = 1.1782595E-02;   Omega22(3,1,2) = -2.1246301E-01;  Omega22(3,1,3) = 8.4561598E+01;
    Omega22(3,2,0) = -1.4719259E-03;  Omega22(3,2,1) = 1.8446968E-02;   Omega22(3,2,2) = -2.6460411E-01;  Omega22(3,2,3) = 1.0911124E+02;
    Omega22(3,3,0) = -7.7439615E-03;  Omega22(3,3,1) = 1.7129007E-01;   Omega22(3,3,2) = -1.4809088E+00;  Omega22(3,3,3) = 2.1284951E+03;
    Omega22(3,4,0) = -5.0478143E-03;  Omega22(3,4,1) = 1.0236186E-01;   Omega22(3,4,2) = -9.0058935E-01;  Omega22(3,4,3) = 4.4472565E+02;
    //O
    Omega22(4,0,0) = -8.3110691E-03;  Omega22(4,0,1) = 1.9617877E-01;   Omega22(4,0,2) = -1.7205427E+00;  Omega22(4,0,3) = 4.0812829E+03;
    Omega22(4,1,0) = -3.7969686E-03;  Omega22(4,1,1) = 7.6789981E-02;   Omega22(4,1,2) = -7.3056809E-01;  Omega22(4,1,3) = 3.3958171E+02;
    Omega22(4,2,0) = -1.0066279E-03;  Omega22(4,2,1) = 1.1029264E-02;   Omega22(4,2,2) = -2.0671266E-01;  Omega22(4,2,3) = 8.2644384E+01;
    Omega22(4,3,0) = -5.0478143E-03;  Omega22(4,3,1) = 1.0236186E-01;   Omega22(4,3,2) = -9.0058935E-01;  Omega22(4,3,3) = 4.4472565E+02;
    Omega22(4,4,0) = -4.2451096E-03;  Omega22(4,4,1) = 9.6820337E-02;   Omega22(4,4,2) = -9.9770795E-01;  Omega22(4,4,3) = 8.3320644E+02;

    // Creation/Destruction (+1/-1), Index of monoatomic reactants
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  1; CatRecombTable(0,1) = 3;
    CatRecombTable(1,0) =  1; CatRecombTable(1,1) = 4;
    CatRecombTable(2,0) =  0; CatRecombTable(2,1) = 0;
    CatRecombTable(3,0) = -1; CatRecombTable(3,1) = 3;
    CatRecombTable(4,0) = -1; CatRecombTable(4,1) = 4;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0241;
      mu_ref[0] = 1.716E-5;
      Sm_ref[0] = 111.0;
      Sk_ref[0] = 194.0;
    }

  } else if (gas_model == "AIR-7"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 7) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }

    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 22;
    ionization  = true;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 0.0; // e-
    RotationModes[1] = 2.0; // N2
    RotationModes[2] = 2.0; // O2
    RotationModes[3] = 2.0; // NO
    RotationModes[4] = 0.0; // N
    RotationModes[5] = 0.0; // O
    RotationModes[6] = 2.0; // NO+

    // Molar mass [kg/kmol]
    MolarMass[0] = 5.4858E-04;      // e-
    MolarMass[1] = 2.0*14.0067;     // N2
    MolarMass[2] = 2.0*15.9994;     // O2
    MolarMass[3] = 14.0067+15.9994; // NO
    MolarMass[4] = 14.0067;         // N
    MolarMass[5] = 15.9994;         // O
    MolarMass[6] = 14.0067+15.9994; // NO+

    //Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;    // e-
    CharVibTemp[1] = 3395.0; // N2
    CharVibTemp[2] = 2239.0; // O2
    CharVibTemp[3] = 2817.0; // NO
    CharVibTemp[4] = 0.0;    // N
    CharVibTemp[5] = 0.0;    // O
    CharVibTemp[6] = 2817.0; // NO+

    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;    // e-
    Enthalpy_Formation[1] = 0.0;    // N2
    Enthalpy_Formation[2] = 0.0;    // O2
    Enthalpy_Formation[3] = 3.0E6;  // NO
    Enthalpy_Formation[4] = 3.36E7; // N
    Enthalpy_Formation[5] = 1.54E7; // O
    Enthalpy_Formation[6] = 3.28E7; // NO+

    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    Ref_Temperature[5] = 0.0;
    Ref_Temperature[6] = 0.0;

    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 0.00E+0;   Blottner(0,1) =  0.00E+0;  Blottner(0,2) = -1.20E1;  // e-
    Blottner(1,0) = 2.68E-2;   Blottner(1,1) =  3.18E-1;  Blottner(1,2) = -1.13E1;  // N2
    Blottner(2,0) = 4.49E-2;   Blottner(2,1) = -8.26E-2;  Blottner(2,2) = -9.20E0;  // O2
    Blottner(3,0) = 4.36E-2;   Blottner(3,1) = -3.36E-2;  Blottner(3,2) = -9.58E0;  // NO
    Blottner(4,0) = 1.16E-2;   Blottner(4,1) =  6.03E-1;  Blottner(4,2) = -1.24E1;  // N
    Blottner(5,0) = 2.03E-2;   Blottner(5,1) =  4.29E-1;  Blottner(5,2) = -1.16E1;  // O
    Blottner(6,0) = 3.02E-1;   Blottner(6,1) =  -3.50E0;  Blottner(6,2) = -3.74E0;  // NO+

    // Number of electron states
    nElStates[0] = 1;  // e-
    nElStates[1] = 15; // N2
    nElStates[2] = 7;  // O2
    nElStates[3] = 16; // NO
    nElStates[4] = 3;  // N
    nElStates[5] = 5;  // O
    nElStates[6] = 8;  // NO+

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    // e: 1 state
    CharElTemp(0,0) = 0.000000000000000E+00;
    ElDegeneracy(0,0) = 1;

    //N2: 15 states
    CharElTemp(1,0)  = 0.000000000000000E+00;
    CharElTemp(1,1)  = 7.223156514095200E+04;
    CharElTemp(1,2)  = 8.577862640384000E+04;
    CharElTemp(1,3)  = 8.605026716160000E+04;
    CharElTemp(1,4)  = 9.535118627874400E+04;
    CharElTemp(1,5)  = 9.805635702203200E+04;
    CharElTemp(1,6)  = 9.968267656935200E+04;
    CharElTemp(1,7)  = 1.048976467715200E+05;
    CharElTemp(1,8)  = 1.116489555200000E+05;
    CharElTemp(1,9)  = 1.225836470400000E+05;
    CharElTemp(1,10) = 1.248856873600000E+05;
    CharElTemp(1,11) = 1.282476158188320E+05;
    CharElTemp(1,12) = 1.338060936000000E+05;
    CharElTemp(1,13) = 1.404296391107200E+05;
    CharElTemp(1,14) = 1.504958859200000E+05;
    ElDegeneracy(1,0)  = 1;
    ElDegeneracy(1,1)  = 3;
    ElDegeneracy(1,2)  = 6;
    ElDegeneracy(1,3)  = 6;
    ElDegeneracy(1,4)  = 3;
    ElDegeneracy(1,5)  = 1;
    ElDegeneracy(1,6)  = 2;
    ElDegeneracy(1,7)  = 2;
    ElDegeneracy(1,8)  = 5;
    ElDegeneracy(1,9)  = 1;
    ElDegeneracy(1,10) = 6;
    ElDegeneracy(1,11) = 6;
    ElDegeneracy(1,12) = 10;
    ElDegeneracy(1,13) = 6;
    ElDegeneracy(1,14) = 6;
    // O2: 7 states
    CharElTemp(2,0) = 0.000000000000000E+00;
    CharElTemp(2,1) = 1.139156019700800E+04;
    CharElTemp(2,2) = 1.898473947826400E+04;
    CharElTemp(2,3) = 4.755973576639200E+04;
    CharElTemp(2,4) = 4.991242097343200E+04;
    CharElTemp(2,5) = 5.092268575561600E+04;
    CharElTemp(2,6) = 7.189863255967200E+04;
    ElDegeneracy(2,0) = 3;
    ElDegeneracy(2,1) = 2;
    ElDegeneracy(2,2) = 1;
    ElDegeneracy(2,3) = 1;
    ElDegeneracy(2,4) = 6;
    ElDegeneracy(2,5) = 3;
    ElDegeneracy(2,6) = 3;
    // NO: 16 states
    CharElTemp(3,0)  = 0.000000000000000E+00;
    CharElTemp(3,1)  = 5.467345760000000E+04;
    CharElTemp(3,2)  = 6.317139627802400E+04;
    CharElTemp(3,3)  = 6.599450342445600E+04;
    CharElTemp(3,4)  = 6.906120960000000E+04;
    CharElTemp(3,5)  = 7.049998480000000E+04;
    CharElTemp(3,6)  = 7.491055017560000E+04;
    CharElTemp(3,7)  = 7.628875293968000E+04;
    CharElTemp(3,8)  = 8.676188537552000E+04;
    CharElTemp(3,9)  = 8.714431182368000E+04;
    CharElTemp(3,10) = 8.886077063728000E+04;
    CharElTemp(3,11) = 8.981755614528000E+04;
    CharElTemp(3,12) = 8.988445919208000E+04;
    CharElTemp(3,13) = 9.042702132000000E+04;
    CharElTemp(3,14) = 9.064283760000000E+04;
    CharElTemp(3,15) = 9.111763341600000E+04;
    ElDegeneracy(3,0)  = 4;
    ElDegeneracy(3,1)  = 8;
    ElDegeneracy(3,2)  = 2;
    ElDegeneracy(3,3)  = 4;
    ElDegeneracy(3,4)  = 4;
    ElDegeneracy(3,5)  = 4;
    ElDegeneracy(3,6)  = 4;
    ElDegeneracy(3,7)  = 2;
    ElDegeneracy(3,8)  = 4;
    ElDegeneracy(3,9)  = 2;
    ElDegeneracy(3,10) = 4;
    ElDegeneracy(3,11) = 4;
    ElDegeneracy(3,12) = 2;
    ElDegeneracy(3,13) = 2;
    ElDegeneracy(3,14) = 2;
    ElDegeneracy(3,15) = 4;
    // N: 3 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.766469645581980E+04;
    CharElTemp(4,2) = 4.149309313560210E+04;
    ElDegeneracy(4,0)= 4;
    ElDegeneracy(4,1)= 10;
    ElDegeneracy(4,2)= 6;
    // O: 5 states
    CharElTemp(5,0) = 0.000000000000000E+00;
    CharElTemp(5,1) = 2.277077570280000E+02;
    CharElTemp(5,2) = 3.265688785704000E+02;
    CharElTemp(5,3) = 2.283028632262240E+04;
    CharElTemp(5,4) = 4.861993036434160E+04;
    ElDegeneracy(5,0) = 5;
    ElDegeneracy(5,1) = 3;
    ElDegeneracy(5,2) = 1;
    ElDegeneracy(5,3) = 5;
    ElDegeneracy(5,4) = 1;
    // NO+: 8 states
    CharElTemp(6,0) = 0.000000000000000E+00;
    CharElTemp(6,1) = 7.508967768800000E+04;
    CharElTemp(6,2) = 8.525462447600000E+04;
    CharElTemp(6,3) = 8.903572570160000E+04;
    CharElTemp(6,4) = 9.746982592400000E+04;
    CharElTemp(6,5) = 1.000553049584000E+05;
    CharElTemp(6,6) = 1.028033655904000E+05;
    CharElTemp(6,7) = 1.057138639424800E+05;
    ElDegeneracy(6,0) = 1;
    ElDegeneracy(6,1) = 3;
    ElDegeneracy(6,2) = 6;
    ElDegeneracy(6,3) = 6;
    ElDegeneracy(6,4) = 3;
    ElDegeneracy(6,5) = 1;
    ElDegeneracy(6,6) = 2;
    ElDegeneracy(6,7) = 2;

    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions(0,0,0)=1;    Reactions(0,0,1)=1;   Reactions(0,0,2)=nSpecies;    Reactions(0,1,0)=4;   Reactions(0,1,1)=4;   Reactions(0,1,2) =1;
    Reactions(1,0,0)=1;    Reactions(1,0,1)=2;   Reactions(1,0,2)=nSpecies;    Reactions(1,1,0)=4;   Reactions(1,1,1)=4;   Reactions(1,1,2) =2;
    Reactions(2,0,0)=1;    Reactions(2,0,1)=3;   Reactions(2,0,2)=nSpecies;    Reactions(2,1,0)=4;   Reactions(2,1,1)=4;   Reactions(2,1,2) =3;
    Reactions(3,0,0)=1;    Reactions(3,0,1)=4;   Reactions(3,0,2)=nSpecies;    Reactions(3,1,0)=4;   Reactions(3,1,1)=4;   Reactions(3,1,2) =4;
    Reactions(4,0,0)=1;    Reactions(4,0,1)=5;   Reactions(4,0,2)=nSpecies;    Reactions(4,1,0)=4;   Reactions(4,1,1)=4;   Reactions(4,1,2) =5;
    Reactions(5,0,0)=1;    Reactions(5,0,1)=6;   Reactions(5,0,2)=nSpecies;    Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2) =6;
    // O2 dissociation
    Reactions(6,0,0)=2;    Reactions(6,0,1)=1;   Reactions(6,0,2)=nSpecies;    Reactions(6,1,0)=5;   Reactions(6,1,1)=5;   Reactions(6,1,2) =1;
    Reactions(7,0,0)=2;    Reactions(7,0,1)=2;   Reactions(7,0,2)=nSpecies;    Reactions(7,1,0)=5;   Reactions(7,1,1)=5;   Reactions(7,1,2) =2;
    Reactions(8,0,0)=2;    Reactions(8,0,1)=3;   Reactions(8,0,2)=nSpecies;    Reactions(8,1,0)=5;   Reactions(8,1,1)=5;   Reactions(8,1,2) =3;
    Reactions(9,0,0)=2;    Reactions(9,0,1)=4;   Reactions(9,0,2)=nSpecies;    Reactions(9,1,0)=5;   Reactions(9,1,1)=5;   Reactions(9,1,2) =4;
    Reactions(10,0,0)=2;   Reactions(10,0,1)=5;  Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=5;  Reactions(10,1,1)=5;  Reactions(10,1,2) =5;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=6;  Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=5;  Reactions(11,1,1)=5;  Reactions(11,1,2) =6;
    // NO dissociation
    Reactions(12,0,0)=3;   Reactions(12,0,1)=1;  Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=4;  Reactions(12,1,1)=5;  Reactions(12,1,2) =1;
    Reactions(13,0,0)=3;   Reactions(13,0,1)=2;  Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=4;  Reactions(13,1,1)=5;  Reactions(13,1,2) =2;
    Reactions(14,0,0)=3;   Reactions(14,0,1)=3;  Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=4;  Reactions(14,1,1)=5;  Reactions(14,1,2) =3;
    Reactions(15,0,0)=3;   Reactions(15,0,1)=4;  Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=4;  Reactions(15,1,1)=5;  Reactions(15,1,2) =4;
    Reactions(16,0,0)=3;   Reactions(16,0,1)=5;  Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=4;  Reactions(16,1,1)=5;  Reactions(16,1,2) =5;
    Reactions(17,0,0)=3;   Reactions(17,0,1)=6;  Reactions(17,0,2)=nSpecies;   Reactions(17,1,0)=4;  Reactions(17,1,1)=5;  Reactions(17,1,2) =6;
    // N2 + O -> NO + N
    Reactions(18,0,0)=1;   Reactions(18,0,1)=5;  Reactions(18,0,2)=nSpecies;   Reactions(18,1,0)=3;  Reactions(18,1,1)=4;  Reactions(18,1,2)= nSpecies;
    // NO + O -> O2 + N
    Reactions(19,0,0)=3;   Reactions(19,0,1)=5;  Reactions(19,0,2)=nSpecies;   Reactions(19,1,0)=2;  Reactions(19,1,1)=4;  Reactions(19,1,2)= nSpecies;
    //N + O -> NO+ + e
    Reactions(20,0,0)=4;   Reactions(20,0,1)=5;  Reactions(20,0,2)=nSpecies;   Reactions(20,1,0)=6;  Reactions(20,1,1)=0;  Reactions(20,1,2)= nSpecies;
    //N2 + e -> N + N + e
    Reactions(21,0,0)=1;   Reactions(21,0,1)=0;  Reactions(21,0,2)=nSpecies;   Reactions(21,1,0)=4;  Reactions(21,1,1)=4;  Reactions(21,1,2)= 0;

    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 7.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 2.0E21;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 1.0E22;
    ArrheniusCoefficient[11] = 2.0E21;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 5.0E15;
    ArrheniusCoefficient[14] = 5.0E15;
    ArrheniusCoefficient[15] = 1.1E17;
    ArrheniusCoefficient[16] = 1.1E17;
    ArrheniusCoefficient[17] = 5.0E15;
    ArrheniusCoefficient[18] = 6.4E17;
    ArrheniusCoefficient[19] = 8.4E12;
    ArrheniusCoefficient[20] = 5.3E12;
    ArrheniusCoefficient[21] = 3.0E24;

    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.60;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = -1.50;
    ArrheniusEta[11] = -1.50;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = 0.0;
    ArrheniusEta[16] = 0.0;
    ArrheniusEta[17] = 0.0;
    ArrheniusEta[18] = -1.0;
    ArrheniusEta[19] = 0.0;
    ArrheniusEta[20] = 0.0;
    ArrheniusEta[21] = -1.60;

    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 113200.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10]  = 59500.0;
    ArrheniusTheta[11]  = 59500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 75500.0;
    ArrheniusTheta[16] = 75500.0;
    ArrheniusTheta[17] = 75500.0;
    ArrheniusTheta[18] = 38400.0;
    ArrheniusTheta[19] = 19450.0;
    ArrheniusTheta[20] = 31900.0;
    ArrheniusTheta[21] = 113200.0;

    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:         a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:             a = 1,   b = 0
    //   Associative ion...    a = 1,   b = 0  ???
    //   E Impact dissociation a = 0,   b = 1
    //   E Impact ionization:  a = 0,   b = 1
    //
    // Backward Reactions
    //   Dissociation:           a = 1,   b = 0
    //   Exchange:               a = 1,   b = 0
    //   Associative  ion...     a = 0.5, b = 0.5
    //   E Impact ionization:    a = 0,   b = 1
    //   E Impact dissocitation: a = 0.5, b = 0.5 ???
    //   N2 impact dissociation: a = 0,   b = 1
    //   Others:                 a = 1,   b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;   Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;   Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;   Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;   Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;   Tcb_b[4] = 0;
    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;   Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;   Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;   Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;   Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;   Tcb_b[9] = 0;
    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;   Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;   Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;   Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;   Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;   Tcb_b[14] = 0;
    Tcf_a[15] = 0.5; Tcf_b[15] = 0.5; Tcb_a[15] = 1;   Tcb_b[15] = 0;
    Tcf_a[16] = 0.5; Tcf_b[16] = 0.5; Tcb_a[16] = 1;   Tcb_b[16] = 0;
    Tcf_a[17] = 0.5; Tcf_b[17] = 0.5; Tcb_a[17] = 1;   Tcb_b[17] = 0;
    Tcf_a[18] = 1.0; Tcf_b[18] = 0.0; Tcb_a[18] = 1;   Tcb_b[18] = 0;
    Tcf_a[19] = 1.0; Tcf_b[19] = 0.0; Tcb_a[19] = 1;   Tcb_b[19] = 0;
    Tcf_a[20] = 1.0; Tcf_b[20] = 0.0; Tcb_a[20] = 0.5; Tcb_b[20] = 0.5;
    Tcf_a[21] = 0.0; Tcf_b[21] = 1.0; Tcb_a[21] = 0;   Tcb_b[21] = 1;

    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3

    // Omega^(1,1) ----------------------
    Omega11(0,0,0) = -1.000000E+00;  Omega11(0,0,1) = -1.000000E+00;  Omega11(0,0,2) = -1.000000E+00;  Omega11(0,0,3) = -1.000000E+00;
    Omega11(0,1,0) = -1.0525124E-02; Omega11(0,1,1) = 1.3498950E-01;  Omega11(0,1,2) = 1.2524805E-01;  Omega11(0,1,3) = 1.5066506E-01;
    Omega11(0,2,0) = 2.3527001E-02;  Omega11(0,2,1) = -6.9632323E-01; Omega11(0,2,2) = 6.8035475E+00;  Omega11(0,2,3) = 1.8335509E-09;
    Omega11(0,3,0) = 1.0414818E-01;  Omega11(0,3,1) = -2.8369126E+00; Omega11(0,3,2) = 2.5323135E+01;  Omega11(0,3,3) = 7.7138358E-32;
    Omega11(0,4,0) = 0.0000000E+00;  Omega11(0,4,1) = 1.6554247E-01;  Omega11(0,4,2) = -3.4986344E+00; Omega11(0,4,3) = 5.9268038E+08;
    Omega11(0,5,0) = 9.9865506E-03;  Omega11(0,5,1) = -2.7407431E-01; Omega11(0,5,2) = 2.6561032E+00;  Omega11(0,5,3) = 4.3080676E-04;
    Omega11(0,6,0) = 1.0000000E+00;  Omega11(0,6,1) = 1.0000000E+00;  Omega11(0,6,2) = 1.0000000E+00;  Omega11(0,6,3) = 1.0000000E+00;
    //N2
    Omega11(1,0,0) = -1.0525124E-02; Omega11(1,0,1) = 1.3498950E-01;  Omega11(1,0,2) = 1.2524805E-01;  Omega11(1,0,3) = 1.5066506E-01;
    Omega11(1,1,0) = -6.0614558E-03; Omega11(1,1,1) = 1.2689102E-01;  Omega11(1,1,2) = -1.0616948E+00; Omega11(1,1,3) = 8.0955466E+02;
    Omega11(1,2,0) = -3.7959091E-03; Omega11(1,2,1) = 9.5708295E-02;  Omega11(1,2,2) = -1.0070611E+00; Omega11(1,2,3) = 8.9392313E+02;
    Omega11(1,3,0) = -1.9295666E-03; Omega11(1,3,1) = 2.7995735E-02;  Omega11(1,3,2) = -3.1588514E-01; Omega11(1,3,3) = 1.2880734E+02;
    Omega11(1,4,0) = -1.0796249E-02; Omega11(1,4,1) = 2.2656509E-01;  Omega11(1,4,2) = -1.7910602E+00; Omega11(1,4,3) = 4.0455218E+03;
    Omega11(1,5,0) = -2.7244269E-03; Omega11(1,5,1) = 6.9587171E-02;  Omega11(1,5,2) = -7.9538667E-01; Omega11(1,5,3) = 4.0673730E+02;
    Omega11(1,6,0) = 0.0000000E+00;  Omega11(1,6,1) = 9.1205839E-02;  Omega11(1,6,2) = -1.8728231E+00; Omega11(1,6,3) = 2.4432020E+05;
    //O2
    Omega11(2,0,0) = 2.3527001E-02;  Omega11(2,0,1) = -6.9632323E-01; Omega11(2,0,2) = 6.8035475E+00;  Omega11(2,0,3) = 1.8335509E-09;
    Omega11(2,1,0) = -3.7959091E-03; Omega11(2,1,1) = 9.5708295E-02;  Omega11(2,1,2) = -1.0070611E+00; Omega11(2,1,3) = 8.9392313E+02;
    Omega11(2,2,0) = -8.0682650E-04; Omega11(2,2,1) = 1.6602480E-02;  Omega11(2,2,2) = -3.1472774E-01; Omega11(2,2,3) = 1.4116458E+02;
    Omega11(2,3,0) = -6.4433840E-04; Omega11(2,3,1) = 8.5378580E-03;  Omega11(2,3,2) = -2.3225102E-01; Omega11(2,3,3) = 1.1371608E+02;
    Omega11(2,4,0) = -1.1453028E-03; Omega11(2,4,1) = 1.2654140E-02;  Omega11(2,4,2) = -2.2435218E-01; Omega11(2,4,3) = 7.7201588E+01;
    Omega11(2,5,0) = -4.8405803E-03; Omega11(2,5,1) = 1.0297688E-01;  Omega11(2,5,2) = -9.6876576E-01; Omega11(2,5,3) = 6.1629812E+02;
    Omega11(2,6,0) = -3.7822765E-03; Omega11(2,6,1) = 1.7967016E-01;  Omega11(2,6,2) = -2.5409098E+00; Omega11(2,6,3) = 1.1840435E+06;
    //NO
    Omega11(3,0,0) = 1.0414818E-01;  Omega11(3,0,1) = -2.8369126E+00; Omega11(3,0,2) = 2.5323135E+01;  Omega11(3,0,3) = 7.7138358E-32;
    Omega11(3,1,0) = -1.9295666E-03; Omega11(3,1,1) = 2.7995735E-02;  Omega11(3,1,2) = -3.1588514E-01; Omega11(3,1,3) = 1.2880734E+02;
    Omega11(3,2,0) = -6.4433840E-04; Omega11(3,2,1) = 8.5378580E-03;  Omega11(3,2,2) = -2.3225102E-01; Omega11(3,2,3) = 1.1371608E+02;
    Omega11(3,3,0) = -0.0000000E+00; Omega11(3,3,1) = -1.1056066E-02; Omega11(3,3,2) = -5.9216250E-02; Omega11(3,3,3) = 7.2542367E+01;
    Omega11(3,4,0) = -1.5770918E-03; Omega11(3,4,1) = 1.9578381E-02;  Omega11(3,4,2) = -2.7873624E-01; Omega11(3,4,3) = 9.9547944E+01;
    Omega11(3,5,0) = -1.0885815E-03; Omega11(3,5,1) = 1.1883688E-02;  Omega11(3,5,2) = -2.1844909E-01; Omega11(3,5,3) = 7.5512560E+01;
    Omega11(3,6,0) = -8.1158474E-03; Omega11(3,6,1) = 2.1474280E-01;  Omega11(3,6,2) = -2.0148450E+00; Omega11(3,6,3) = 6.2986385E+04;
    //N
    Omega11(4,0,0) = 0.0000000E+00;  Omega11(4,0,1) = 1.6554247E-01;  Omega11(4,0,2) = -3.4986344E+00; Omega11(4,0,3) = 5.9268038E+08;
    Omega11(4,1,0) = -1.0796249E-02; Omega11(4,1,1) = 2.2656509E-01;  Omega11(4,1,2) = -1.7910602E+00; Omega11(4,1,3) = 4.0455218E+03;
    Omega11(4,2,0) = -1.1453028E-03; Omega11(4,2,1) = 1.2654140E-02;  Omega11(4,2,2) = -2.2435218E-01; Omega11(4,2,3) = 7.7201588E+01;
    Omega11(4,3,0) = -1.5770918E-03; Omega11(4,3,1) = 1.9578381E-02;  Omega11(4,3,2) = -2.7873624E-01; Omega11(4,3,3) = 9.9547944E+01;
    Omega11(4,4,0) = -9.6083779E-03; Omega11(4,4,1) = 2.0938971E-01;  Omega11(4,4,2) = -1.7386904E+00; Omega11(4,4,3) = 3.3587983E+03;
    Omega11(4,5,0) = -7.8147689E-03; Omega11(4,5,1) = 1.6792705E-01;  Omega11(4,5,2) = -1.4308628E+00; Omega11(4,5,3) = 1.6628859E+03;
    Omega11(4,6,0) = -1.9605234E-02; Omega11(4,6,1) = 5.5570872E-01;  Omega11(4,6,2) = -5.4285702E+00; Omega11(4,6,3) = 1.3574446E+09;
    //O
    Omega11(5,0,0) = 9.9865506E-03;  Omega11(5,0,1) = -2.7407431E-01; Omega11(5,0,2) = 2.6561032E+00;  Omega11(5,0,3) = 4.3080676E-04;
    Omega11(5,1,0) = -2.7244269E-03; Omega11(5,1,1) = 6.9587171E-02;  Omega11(5,1,2) = -7.9538667E-01; Omega11(5,1,3) = 4.0673730E+02;
    Omega11(5,2,0) = -4.8405803E-03; Omega11(5,2,1) = 1.0297688E-01;  Omega11(5,2,2) = -9.6876576E-01; Omega11(5,2,3) = 6.1629812E+02;
    Omega11(5,3,0) = -1.0885815E-03; Omega11(5,3,1) = 1.1883688E-02;  Omega11(5,3,2) = -2.1844909E-01; Omega11(5,3,3) = 7.5512560E+01;
    Omega11(5,4,0) = -7.8147689E-03; Omega11(5,4,1) = 1.6792705E-01;  Omega11(5,4,2) = -1.4308628E+00; Omega11(5,4,3) = 1.6628859E+03;
    Omega11(5,5,0) = -6.4040535E-03; Omega11(5,5,1) = 1.4629949E-01;  Omega11(5,5,2) = -1.3892121E+00; Omega11(5,5,3) = 2.0903441E+03;
    Omega11(5,6,0) = -1.6409054E-02; Omega11(5,6,1) = 4.6352852E-01;  Omega11(5,6,2) = -4.5479735E+00; Omega11(5,6,3) = 7.4250671E+07;
    //NO+
    Omega11(6,0,0) = 1.0000000E+00;  Omega11(6,0,1) = 1.0000000E+00;  Omega11(6,0,2) = 1.0000000E+00;  Omega11(6,0,3) = 1.0000000E+00;
    Omega11(6,1,0) = 0.0000000E+00;  Omega11(6,1,1) = 9.1205839E-02;  Omega11(6,1,2) = -1.8728231E+00; Omega11(6,1,3) = 2.4432020E+05;
    Omega11(6,2,0) = -3.7822765E-03; Omega11(6,2,1) = 1.7967016E-01;  Omega11(6,2,2) = -2.5409098E+00; Omega11(6,2,3) = 1.1840435E+06;
    Omega11(6,3,0) = -8.1158474E-03; Omega11(6,3,1) = 2.1474280E-01;  Omega11(6,3,2) = -2.0148450E+00; Omega11(6,3,3) = 6.2986385E+04;
    Omega11(6,4,0) = -1.9605234E-02; Omega11(6,4,1) = 5.5570872E-01;  Omega11(6,4,2) = -5.4285702E+00; Omega11(6,4,3) = 1.3574446E+09;
    Omega11(6,5,0) = -1.6409054E-02; Omega11(6,5,1) = 4.6352852E-01;  Omega11(6,5,2) = -4.5479735E+00; Omega11(6,5,3) = 7.4250671E+07;
    Omega11(6,6,0) = -1.000000E+00;  Omega11(6,6,1) = -1.000000E+00;  Omega11(6,6,2) = -1.000000E+00;  Omega11(6,6,3) = -1.000000E+00;

    // Omega^(2,2) ----------------------
    Omega22(0,0,0) = -1.000000E+00;  Omega22(0,0,1) = -1.000000E+00;  Omega22(0,0,2) = -1.000000E+00;  Omega22(0,0,3) = -1.000000E+00;
    Omega22(0,1,0) = -4.2254948E-03; Omega22(0,1,1) = -5.2965163E-02; Omega22(0,1,2) = 1.9157708E+00;  Omega22(0,1,3) = 6.3263309E-04;
    Omega22(0,2,0) = 9.6744867E-03;  Omega22(0,2,1) = -3.3759583E-01; Omega22(0,2,2) = 3.7952121E+00;  Omega22(0,2,3) = 6.8468036E-06;
    Omega22(0,3,0) = 0.0000000E+00;  Omega22(0,3,1) = 5.4444485E-02;  Omega22(0,3,2) = -1.2854128E+00; Omega22(0,3,3) = 1.3857556E+04;
    Omega22(0,4,0) = -1.0903638E-01; Omega22(0,4,1) = 2.8678381E+00;  Omega22(0,4,2) = -2.5297550E+01; Omega22(0,4,3) = 3.4838798E+33;
    Omega22(0,5,0) = -1.7924100E-02; Omega22(0,5,1) = 4.0402656E-01;  Omega22(0,5,2) = -2.6712374E+00; Omega22(0,5,3) = 4.1447669E+02;
    Omega22(0,6,0) = 1.0000000E+00;  Omega22(0,6,1) = 1.0000000E+00;  Omega22(0,6,2) = 1.0000000E+00;  Omega22(0,6,3) = 1.0000000E+00;
    //N2
    Omega22(1,0,0) = -4.2254948E-03; Omega22(1,0,1) = -5.2965163E-02; Omega22(1,0,2) = 1.9157708E+00;  Omega22(1,0,3) = 6.3263309E-04;
    Omega22(1,1,0) = -7.6303990E-03; Omega22(1,1,1) = 1.6878089E-01;  Omega22(1,1,2) = -1.4004234E+00; Omega22(1,1,3) = 2.1427708E+03;
    Omega22(1,2,0) = -8.0457321E-03; Omega22(1,2,1) = 1.9228905E-01;  Omega22(1,2,2) = -1.7102854E+00; Omega22(1,2,3) = 5.2213857E+03;
    Omega22(1,3,0) = -6.8237776E-03; Omega22(1,3,1) = 1.4360616E-01;  Omega22(1,3,2) = -1.1922240E+00; Omega22(1,3,3) = 1.2433086E+03;
    Omega22(1,4,0) = -8.3493693E-03; Omega22(1,4,1) = 1.7808911E-01;  Omega22(1,4,2) = -1.4466155E+00; Omega22(1,4,3) = 1.9324210E+03;
    Omega22(1,5,0) = -8.3110691E-03; Omega22(1,5,1) = 1.9617877E-01;  Omega22(1,5,2) = -1.7205427E+00; Omega22(1,5,3) = 4.0812829E+03;
    Omega22(1,6,0) = 0.0000000E+00;  Omega22(1,6,1) = 8.5112236E-02;  Omega22(1,6,2) = -1.7460044E+00; Omega22(1,6,3) = 1.4498969E+05;
    //O2
    Omega22(2,0,0) = 9.6744867E-03;  Omega22(2,0,1) = -3.3759583E-01; Omega22(2,0,2) = 3.7952121E+00;  Omega22(2,0,3) = 6.8468036E-06;
    Omega22(2,1,0) = -8.0457321E-03; Omega22(2,1,1) = 1.9228905E-01;  Omega22(2,1,2) = -1.7102854E+00; Omega22(2,1,3) = 5.2213857E+03;
    Omega22(2,2,0) = -6.2931612E-03; Omega22(2,2,1) = 1.4624645E-01;  Omega22(2,2,2) = -1.3006927E+00; Omega22(2,2,3) = 1.8066892E+03;
    Omega22(2,3,0) = -6.8508672E-03; Omega22(2,3,1) = 1.5524564E-01;  Omega22(2,3,2) = -1.3479583E+00; Omega22(2,3,3) = 2.0037890E+03;
    Omega22(2,4,0) = -1.0608832E-03; Omega22(2,4,1) = 1.1782595E-02;  Omega22(2,4,2) = -2.1246301E-01; Omega22(2,4,3) = 8.4561598E+01;
    Omega22(2,5,0) = -3.7969686E-03; Omega22(2,5,1) = 7.6789981E-02;  Omega22(2,5,2) = -7.3056809E-01; Omega22(2,5,3) = 3.3958171E+02;
    Omega22(2,6,0) = 0.0000000E+00;  Omega22(2,6,1) = 8.4737359E-02;  Omega22(2,6,2) = -1.7290488E+00; Omega22(2,6,3) = 1.2485194E+05;
    //NO
    Omega22(3,0,0) = 0.0000000E+00;  Omega22(0,3,1) = 5.4444485E-02;  Omega22(0,3,2) = -1.2854128E+00; Omega22(0,3,3) = 1.3857556E+04;
    Omega22(3,1,0) = -6.8237776E-03; Omega22(3,1,1) = 1.4360616E-01;  Omega22(3,1,2) = -1.1922240E+00; Omega22(3,1,3) = 1.2433086E+03;
    Omega22(3,2,0) = -6.8508672E-03; Omega22(3,2,1) = 1.5524564E-01;  Omega22(3,2,2) = -1.3479583E+00; Omega22(3,2,3) = 2.0037890E+03;
    Omega22(3,3,0) = -7.4942466E-03; Omega22(3,3,1) = 1.6626193E-01;  Omega22(3,3,2) = -1.4107027E+00; Omega22(3,3,3) = 2.3097604E+03;
    Omega22(3,4,0) = -1.4719259E-03; Omega22(3,4,1) = 1.8446968E-02;  Omega22(3,4,2) = -2.6460411E-01; Omega22(3,4,3) = 1.0911124E+02;
    Omega22(3,5,0) = -1.0066279E-03; Omega22(3,5,1) = 1.1029264E-02;  Omega22(3,5,2) = -2.0671266E-01; Omega22(3,5,3) = 8.2644384E+01;
    Omega22(3,6,0) = 1.1055777E-02;  Omega22(3,6,1) = -1.6621846E-01; Omega22(3,6,2) = 1.4372166E-01;  Omega22(3,6,3) = 1.3182061E+03;
    //N
    Omega22(4,0,0) = -1.0903638E-01; Omega22(4,0,1) = 2.8678381E+00;  Omega22(4,0,2) = -2.5297550E+01; Omega22(4,0,3) = 3.4838798E+33;
    Omega22(4,1,0) = -8.3493693E-03; Omega22(4,1,1) = 1.7808911E-01;  Omega22(4,1,2) = -1.4466155E+00; Omega22(4,1,3) = 1.9324210E+03;
    Omega22(4,2,0) = -1.0608832E-03; Omega22(4,2,1) = 1.1782595E-02;  Omega22(4,2,2) = -2.1246301E-01; Omega22(4,2,3) = 8.4561598E+01;
    Omega22(4,3,0) = -1.4719259E-03; Omega22(4,3,1) = 1.8446968E-02;  Omega22(4,3,2) = -2.6460411E-01; Omega22(4,3,3) = 1.0911124E+02;
    Omega22(4,4,0) = -7.7439615E-03; Omega22(4,4,1) = 1.7129007E-01;  Omega22(4,4,2) = -1.4809088E+00; Omega22(4,4,3) = 2.1284951E+03;
    Omega22(4,5,0) = -5.0478143E-03; Omega22(4,5,1) = 1.0236186E-01;  Omega22(4,5,2) = -9.0058935E-01; Omega22(4,5,3) = 4.4472565E+02;
    Omega22(4,6,0) = -2.1009546E-02; Omega22(4,6,1) = 5.8910426E-01;  Omega22(4,6,2) = -5.6681361E+00; Omega22(4,6,3) = 2.4486594E+09;
    //O
    Omega22(5,0,0) = -1.7924100E-02; Omega22(5,0,1) = 4.0402656E-01;  Omega22(5,0,2) = -2.6712374E+00; Omega22(5,0,3) = 4.1447669E+02;
    Omega22(5,1,0) = -8.3110691E-03; Omega22(5,1,1) = 1.9617877E-01;  Omega22(5,1,2) = -1.7205427E+00; Omega22(5,1,3) = 4.0812829E+03;
    Omega22(5,2,0) = -3.7969686E-03; Omega22(5,2,1) = 7.6789981E-02;  Omega22(5,2,2) = -7.3056809E-01; Omega22(5,2,3) = 3.3958171E+02;
    Omega22(5,3,0) = -1.0066279E-03; Omega22(5,3,1) = 1.1029264E-02;  Omega22(5,3,2) = -2.0671266E-01; Omega22(5,3,3) = 8.2644384E+01;
    Omega22(5,4,0) = -5.0478143E-03; Omega22(5,4,1) = 1.0236186E-01;  Omega22(5,4,2) = -9.0058935E-01; Omega22(5,4,3) = 4.4472565E+02;
    Omega22(5,5,0) = -4.2451096E-03; Omega22(5,5,1) = 9.6820337E-02;  Omega22(5,5,2) = -9.9770795E-01; Omega22(5,5,3) = 8.3320644E+02;
    Omega22(5,6,0) = -1.5315132E-02; Omega22(5,6,1) = 4.3541627E-01;  Omega22(5,6,2) = -4.2864279E+00; Omega22(5,6,3) = 3.5125207E+07;
    //NO+
    Omega22(6,0,0) = 1.0000000E+00;  Omega22(6,0,1) = 1.0000000E+00;  Omega22(6,0,2) = 1.0000000E+00;  Omega22(6,0,3) = 1.0000000E+00;
    Omega22(6,1,0) = 0.0000000E+00;  Omega22(6,1,1) = 8.5112236E-02;  Omega22(6,1,2) = -1.7460044E+00; Omega22(6,1,3) = 1.4498969E+05;
    Omega22(6,2,0) = 0.0000000E+00;  Omega22(6,2,1) = 8.4737359E-02;  Omega22(6,2,2) = -1.7290488E+00; Omega22(6,2,3) = 1.2485194E+05;
    Omega22(6,3,0) = 1.1055777E-02;  Omega22(6,3,1) = -1.6621846E-01; Omega22(6,3,2) = 1.4372166E-01;  Omega22(6,3,3) = 1.3182061E+03;
    Omega22(6,4,0) = -2.1009546E-02; Omega22(6,4,1) = 5.8910426E-01;  Omega22(6,4,2) = -5.6681361E+00; Omega22(6,4,3) = 2.4486594E+09;
    Omega22(6,5,0) = -1.5315132E-02; Omega22(6,5,1) = 4.3541627E-01;  Omega22(6,5,2) = -4.2864279E+00; Omega22(6,5,3) = 3.5125207E+07;
    Omega22(6,6,0) = -1.000000E+00;  Omega22(6,6,1) = -1.000000E+00;  Omega22(6,6,2) = -1.000000E+00;  Omega22(6,6,3) = -1.000000E+00;

    // Creation/Destruction (+1/-1), Index of monoatomic reactants
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  0; CatRecombTable(0,1) = 1;
    CatRecombTable(1,0) =  1; CatRecombTable(1,1) = 4;
    CatRecombTable(2,0) =  1; CatRecombTable(2,1) = 5;
    CatRecombTable(3,0) =  0; CatRecombTable(3,1) = 1;
    CatRecombTable(4,0) = -1; CatRecombTable(4,1) = 4;
    CatRecombTable(5,0) = -1; CatRecombTable(5,1) = 5;
    CatRecombTable(6,0) =  0; CatRecombTable(6,1) = 1;

    /*--- Values for Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0241;
      mu_ref[0] = 1.716E-5;
      Sm_ref[0] = 111.0;
      Sk_ref[0] = 194.0;
    }
  }

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
}

CSU2TCLib::~CSU2TCLib()= default;

void CSU2TCLib::SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){

  rhos = val_rhos;
  T    = val_temperature;
  Tve  = single_temperature ? val_temperature : val_temperature_ve;

  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density += rhos[iSpecies];

  Pressure = ComputePressure();

}

vector<su2double>& CSU2TCLib::GetSpeciesCvTraRot(){

  if(ionization) Cvtrs[0] = 0.0;

  for (iSpecies = nEl; iSpecies < nHeavy; iSpecies++)
    Cvtrs[iSpecies] = (3.0/2.0 + RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies];

  return Cvtrs;
}

vector<su2double>& CSU2TCLib::ComputeSpeciesCvVibEle(su2double val_T){

  su2double thoTve, exptv, num, num2, num3, denom, Cvvs, Cves;
  unsigned short iElectron = 0;

  /*--- Loop through species ---*/
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- If requesting electron specific heat ---*/
    if (ionization && iSpecies == iElectron) {
      Cvvs = 0.0;
      Cves = 3.0/2.0 * Ru/MolarMass[iSpecies];
    }

    /*--- Heavy particle specific heat ---*/
    else {

      /*--- Vibrational energy ---*/
      if (CharVibTemp[iSpecies] != 0.0) {
        thoTve = CharVibTemp[iSpecies]/val_T;
        exptv = exp(CharVibTemp[iSpecies]/val_T);
        Cvvs  = Ru/MolarMass[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
      } else {
        Cvvs = 0.0;
      }

      /*--- Electronic energy ---*/
      if (nElStates[iSpecies] != 0) {
        num = 0.0; num2 = 0.0;
        denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/val_T);
        num3  = ElDegeneracy[iSpecies][0] * (CharElTemp[iSpecies][0]/(val_T*val_T))*exp(-CharElTemp[iSpecies][0]/val_T);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          thoTve = CharElTemp[iSpecies][iEl]/val_T;
          exptv = exp(-CharElTemp[iSpecies][iEl]/val_T);

          num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exptv;
          denom += ElDegeneracy[iSpecies][iEl] * exptv;
          num2  += ElDegeneracy[iSpecies][iEl] * (thoTve*thoTve) * exptv;
          num3  += ElDegeneracy[iSpecies][iEl] * thoTve/val_T * exptv;
        }
        Cves = Ru/MolarMass[iSpecies] * (num2/denom - num*num3/(denom*denom));
      } else {
        Cves = 0.0;
      }
    }

    Cvves[iSpecies] = Cvvs + Cves;
  }

  return Cvves;

}

vector<su2double>& CSU2TCLib::ComputeMixtureEnergies(){

  su2double Ev, Ee, Ef, num;

  su2double rhoEmix = 0.0;
  su2double rhoEve  = 0.0;
  su2double denom   = 0.0;

  // Electrons
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {

    // Species formation energy
    Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies] * Ref_Temperature[iSpecies];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += rhos[iSpecies]*((3.0/2.0) * Ru/MolarMass[iSpecies] * (Tve - Ref_Temperature[iSpecies]));
  }

  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++){

    // Species formation energy
    Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies];

    // Species vibrational energy
    if (CharVibTemp[iSpecies] != 0.0)
      Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = ElDegeneracy(iSpecies,0) * exp(-CharElTemp(iSpecies,0)/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += ElDegeneracy(iSpecies,iEl) * CharElTemp(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
      denom += ElDegeneracy(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
    }
    Ee = Ru/MolarMass[iSpecies] * (num/denom);

    // Mixture total energy
    rhoEmix += rhos[iSpecies] * ((3.0/2.0+RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies] * (T-Ref_Temperature[iSpecies]) + Ev + Ee + Ef);

    // Mixture vibrational-electronic energy
    rhoEve += rhos[iSpecies] * (Ev + Ee);

  }

  energies[0] = rhoEmix/Density;
  energies[1] = rhoEve/Density;

  return energies;

}

vector<su2double>& CSU2TCLib::ComputeSpeciesEve(su2double val_T, bool vibe_only){

  su2double Ev, Eel, Ef, num, denom;
  unsigned short iElectron = 0;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- Electron species energy ---*/
    if ( ionization && (iSpecies == iElectron)) {

      /*--- Calculate formation energy ---*/
      Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies] * Ref_Temperature[iSpecies];

      /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
      Eel = (3.0/2.0) * Ru/MolarMass[iSpecies] * (val_T - Ref_Temperature[iSpecies]) + Ef;
      Ev  = 0.0;
    }
    /*--- Heavy particle energy ---*/
    else {

      /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
      if (CharVibTemp[iSpecies] != 0.0)
        Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/val_T)-1.0);
      else
        Ev = 0.0;

      /*--- Calculate electronic energy ---*/
      num = 0.0;
      denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/val_T);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
        denom += ElDegeneracy[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
      }
      Eel = Ru/MolarMass[iSpecies] * (num/denom);
    }
    if (vibe_only) {eves[iSpecies] = Ev;}
    else {eves[iSpecies] = Ev + Eel;}
  }

  return eves;
}

vector<su2double>& CSU2TCLib::ComputeNetProductionRates(bool implicit, const su2double *V, const su2double* eve,
                                                        const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                                                        su2double **val_jacobian){

  /*---                          ---*/
  /*--- Nonequilibrium chemistry ---*/
  /*---                          ---*/

  /*--- Initialize variables ---*/
  unsigned short ii, iReaction;
  ws.resize(nSpecies,0.0);
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
    ws[iSpecies] = 0.0;

  /*--- Define artificial chemistry parameters ---*/
  // Note: These parameters artificially increase the rate-controlling reaction
  //       temperature.  This relaxes some of the stiffness in the chemistry
  //       source term.
  const su2double T_min   = 800.0;
  const su2double epsilon = 80;

  /*--- Define preferential dissociation coefficient ---*/
  //alpha = 0.3; //TODO: make this a config option?

  /*--- Loop over all reactions ---*/
  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    /*--- Determine the rate-controlling temperature ---*/
    af = Tcf_a[iReaction];
    bf = Tcf_b[iReaction];
    ab = Tcb_a[iReaction];
    bb = Tcb_b[iReaction];
    if (single_temperature) {
      Trxnf = T;
      Trxnb = T;
    } else {
      Trxnf = pow(T, af)*pow(Tve, bf);
      Trxnb = pow(T, ab)*pow(Tve, bb);
    }

    /*--- Calculate the modified temperature ---*/
    Thf = 0.5 * (Trxnf+T_min + sqrt((Trxnf-T_min)*(Trxnf-T_min)+epsilon*epsilon));
    Thb = 0.5 * (Trxnb+T_min + sqrt((Trxnb-T_min)*(Trxnb-T_min)+epsilon*epsilon));

    /*--- Get the Keq & Arrhenius coefficients ---*/
    ComputeKeqConstants(iReaction);

    /*--- Calculate Keq ---*/
    const su2double Keq = exp(  A[0]*(Thb/1E4) + A[1] + A[2]*log(1E4/Thb)
        + A[3]*(1E4/Thb) + A[4]*(1E4/Thb)*(1E4/Thb) );

    /*--- Calculate rate coefficients ---*/
    kf  = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thf)) * exp(-ArrheniusTheta[iReaction]/Thf);
    kfb = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thb)) * exp(-ArrheniusTheta[iReaction]/Thb);
    kb  = kfb / Keq;

    /*--- Determine production & destruction of each species ---*/
    fwdRxn = 1.0;
    bkwRxn = 1.0;
    for (ii = 0; ii < 3; ii++) {

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if ( iSpecies != nSpecies)
        fwdRxn *= 0.001*rhos[iSpecies]/MolarMass[iSpecies];

      /*--- Products ---*/
      jSpecies = Reactions(iReaction,1,ii);
      if (jSpecies != nSpecies) {
        bkwRxn *= 0.001*rhos[jSpecies]/MolarMass[jSpecies];
      }
    }

    fwdRxn = 1000.0 * kf * fwdRxn;
    bkwRxn = 1000.0 * kb * bkwRxn;

    for (ii = 0; ii < 3; ii++) {

      /*--- Products ---*/
      iSpecies = Reactions(iReaction,1,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] += MolarMass[iSpecies] * (fwdRxn-bkwRxn);

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] -= MolarMass[iSpecies] * (fwdRxn-bkwRxn);
    }

    if (implicit) {
      ChemistryJacobian(iReaction, V, eve, cvve, dTdU, dTvedU, val_jacobian);
    }
  } //iReaction

  return ws;
}

void CSU2TCLib::ChemistryJacobian(unsigned short iReaction, const su2double *V,
                                  const su2double* eve, const su2double *cvve,
                                  const su2double* dTdU, const su2double* dTvedU,
                                  su2double **val_jacobian) {

  unsigned short ii, iVar, jVar, iSpecies;
  unsigned short nEve = nSpecies+nDim+1;
  unsigned short nVar = nSpecies+nDim+2;

  su2double T_min   = 800.0;
  su2double epsilon = 80;

  /*--- Initializing derivative variables ---*/
  dkf.resize(nVar,0.0);      dkb.resize(nVar,0.0);
  dRfok.resize(nVar,0.0);    dRbok.resize(nVar,0.0);
  alphak.resize(nSpecies,0); betak.resize(nSpecies,0);

  for (iVar=0;iVar<nVar;iVar++){
   dkf[iVar]=0.0; dRfok[iVar]=0.0;
   dkb[iVar]=0.0; dRbok[iVar]=0.0;
  }

  for (iSpecies=0;iSpecies<nSpecies;iSpecies++){
   alphak[iSpecies]=0; betak[iSpecies]=0;
  }

  /*--- Extract additional Arrhenius information ---*/
  su2double eta   = ArrheniusEta[iReaction];
  su2double theta = ArrheniusTheta[iReaction];

  /*--- Derivative of modified temperature wrt Trxnf ---*/
  dThf = 0.5 * (1.0 + (Trxnf-T_min)/sqrt((Trxnf-T_min)*(Trxnf-T_min)
                                                  + epsilon*epsilon));
  dThb = 0.5 * (1.0 + (Trxnb-T_min)/sqrt((Trxnb-T_min)*(Trxnb-T_min)
                                                  + epsilon*epsilon));

  /*--- Fwd rate coefficient derivatives ---*/
  coeff = kf * (eta/Thf+theta/(Thf*Thf)) * dThf;
  for (iVar = 0; iVar < nVar; iVar++) {
    if (single_temperature) {
      dkf[iVar] = coeff * dTdU[iVar];
    } else {
      dkf[iVar] = coeff * ( af*Trxnf/T*dTdU[iVar] +
                            bf*Trxnf/Tve*dTvedU[iVar] );
    }
  }

  /*--- Bkwd rate coefficient derivatives ---*/
  coeff = kb * (eta/Thb+theta/(Thb*Thb)) * dThb;
  for (iVar = 0; iVar < nVar; iVar++) {
    if (single_temperature) {
      dkb[iVar] = coeff*dTdU[iVar]
                        - kb*((A[0]*Thb/1E4 - A[2] - A[3]*1E4/Thb
                        - 2*A[4]*(1E4/Thb)*(1E4/Thb))/Thb) * dThb * dTdU[iVar];
    } else {
      dkb[iVar] = coeff*( ab*Trxnb/T*dTdU[iVar] + bb*Trxnb/Tve*dTvedU[iVar])
                        - kb*((A[0]*Thb/1E4 - A[2] - A[3]*1E4/Thb
                        - 2*A[4]*(1E4/Thb)*(1E4/Thb))/Thb) * dThb *
                        ( ab*Trxnb/T*dTdU[iVar] + bb*Trxnb/Tve*dTvedU[iVar]);
    }
  }

  /*--- Rxn rate derivatives ---*/
  for (ii = 0; ii < 3; ii++) {

    /*--- Products ---*/
    iSpecies = Reactions(iReaction,1,ii);
    if (iSpecies != nSpecies)
      betak[iSpecies]++;

    /*--- Reactants ---*/
    iSpecies = Reactions(iReaction,0,ii);
    if (iSpecies != nSpecies)
      alphak[iSpecies]++;
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    // Fwd
    dRfok[iSpecies] =  0.001*alphak[iSpecies]/MolarMass[iSpecies] *
                       pow(0.001*rhos[iSpecies]/MolarMass[iSpecies],
                       max(0, alphak[iSpecies]-1));

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
      if (jSpecies != iSpecies)
        dRfok[iSpecies] *= pow(0.001*rhos[jSpecies]/MolarMass[jSpecies],
                                                       alphak[jSpecies]);
    dRfok[iSpecies] *= 1000.0;

    // Bkw
    dRbok[iSpecies] =  0.001*betak[iSpecies]/MolarMass[iSpecies] *
                       pow(0.001*rhos[iSpecies]/MolarMass[iSpecies],
                       max(0, betak[iSpecies]-1));

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
      if (jSpecies != iSpecies)
        dRbok[iSpecies] *= pow(0.001*rhos[jSpecies]/MolarMass[jSpecies],
                                                        betak[jSpecies]);
    dRbok[iSpecies] *= 1000.0;
  }

  for (ii = 0; ii < 3; ii++) {

    /*--- Products ---*/
    iSpecies = Reactions(iReaction,1,ii);
    if (iSpecies != nSpecies) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[iSpecies][iVar] += MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]); //TODO * Volume;
        if (!single_temperature) {
          val_jacobian[nEve][iVar]   += MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) *
                                                                eve[iSpecies];//TODO * Volume;
        }
      }

      for (jVar = 0; jVar < nVar; jVar++) {
        if (!single_temperature) {
          val_jacobian[nEve][jVar] += MolarMass[iSpecies] * (fwdRxn-bkwRxn)* cvve[iSpecies] *
                                                                               dTvedU[jVar];//TODO * Volume;
        }
      }
    }

    /*--- Reactants ---*/
    iSpecies = Reactions(iReaction,0,ii);
    if (iSpecies != nSpecies) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[iSpecies][iVar] -= MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]);//TODO * Volume;
        if (!single_temperature) {
          val_jacobian[nEve][iVar] -=   MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) *
                                                                                       eve[iSpecies];//TODO * Volume;
        }
      }

      for (jVar = 0; jVar < nVar; jVar++) {
        if (!single_temperature) {
          val_jacobian[nEve][jVar] -= MolarMass[iSpecies] * (fwdRxn-bkwRxn) * cvve[iSpecies] *
                                                                                dTvedU[jVar];//TODO * Volume;
        }
      }
    }
  } // ii
}

void CSU2TCLib::ComputeKeqConstants(unsigned short val_Reaction) {

  unsigned short ii;

  /*--- Acquire database constants from CConfig ---*/
  GetChemistryEquilConstants(val_Reaction);

  /*--- Calculate mixture number density ---*/
  su2double N = 0.0;
  for (iSpecies =0 ; iSpecies < nSpecies; iSpecies++) {
    N += rhos[iSpecies]/MolarMass[iSpecies]*AVOGAD_CONSTANT;
  }

  /*--- Convert number density from 1/m^3 to 1/cm^3 for table look-up ---*/
  N = N*(1E-6);

  /*--- Determine table index based on mixture N ---*/
  unsigned short tbl_offset = 14;
  unsigned short pwr        = floor(log10(N));

  /*--- Bound the interpolation to table limit values ---*/
  unsigned short iIndex = int(pwr) - tbl_offset;
  if (iIndex <= 0) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(0,ii);
    return;
  } if (iIndex >= 5) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(5,ii);
    return;
  }

  /*--- Calculate interpolation denominator terms avoiding pow() ---*/
  su2double tmp1 = 1.0;
  su2double tmp2 = 1.0;
  for (ii = 0; ii < pwr; ii++) {
    tmp1 *= 10.0;
    tmp2 *= 10.0;
  }
  tmp2 *= 10.0;

  /*--- Interpolate ---*/
  for (ii = 0; ii < 5; ii++) {
    A[ii] =  (RxnConstantTable(iIndex+1,ii) - RxnConstantTable(iIndex,ii))
        / (tmp2 - tmp1) * (N - tmp1)
        + RxnConstantTable(iIndex,ii);
  }
}

su2double CSU2TCLib::ComputeEveSourceTerm(){

  if (single_temperature) {
    omega = 0.0;
    return omega;
  }

  /*---                                                                    ---*/
  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  /*---                                                                    ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section

  su2activematrix mu;
  vector<su2double> MolarFrac;

  MolarFrac.resize(nSpecies,0.0);
  mu.resize(nSpecies,nSpecies)=su2double(0.0);

  su2double omegaVT = 0.0;
  su2double omegaCV = 0.0;

  /*--- Calculate mole fractions ---*/
  su2double N    = 0.0;
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc += rhos[iSpecies] / MolarMass[iSpecies];
    N    += rhos[iSpecies] / MolarMass[iSpecies] * AVOGAD_CONSTANT;
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFrac[iSpecies] = (rhos[iSpecies] / MolarMass[iSpecies]) / conc;

  /*--- Compute Eve and Eve* ---*/
  eve_eq = ComputeSpeciesEve(T, true);
  eve    = ComputeSpeciesEve(Tve, true);

  /*--- Loop over species to calculate source term --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Millikan & White relaxation time ---*/
    su2double num   = 0.0;
    su2double denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      mu(iSpecies,jSpecies) = MolarMass[iSpecies]*MolarMass[jSpecies] / (MolarMass[iSpecies] + MolarMass[jSpecies]);
      const su2double A_sr   = 1.16 * 1E-3 * sqrt(mu(iSpecies,jSpecies)) * pow(CharVibTemp[iSpecies], 4.0/3.0);
      const su2double B_sr   = 0.015 * pow(mu(iSpecies,jSpecies), 0.25);
      const su2double tau_sr = 101325.0/Pressure * exp(A_sr*(pow(T,-1.0/3.0) - B_sr) - 18.42);

      num   += MolarFrac[jSpecies];
      denom += MolarFrac[jSpecies] / tau_sr;
    }

    const su2double tauMW = num / denom;

    /*--- Park limiting cross section ---*/
    const su2double Cs    = sqrt((8.0*Ru*T)/(PI_NUMBER*MolarMass[iSpecies]));
    const su2double sig_s = 3E-21*(2.5E9)/(T*T);

    const su2double tauP = 1/(sig_s*Cs*N);

    /*--- Species relaxation time ---*/
    taus[iSpecies] = tauMW + tauP;

    /*--- Add species contribution to residual ---*/
    omegaVT += rhos[iSpecies] * (eve_eq[iSpecies] -
                                 eve[iSpecies]) / taus[iSpecies];
  }

  /*--- Vibrational energy change due to chemical reactions ---*/
  if(!frozen){
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      omegaCV += ws[iSpecies]*eve[iSpecies];
  }

  omega = omegaVT + omegaCV;

  return omega;

}

void CSU2TCLib::GetEveSourceTermJacobian(const su2double *V, const su2double *eve, const su2double *cvve, const su2double *dTdU, const su2double* dTvedU, su2double **val_jacobian){

  if (single_temperature) {
    return;
  }

  unsigned short iVar;
  unsigned short nEv  = nSpecies+nDim+1;
  unsigned short nVar = nSpecies+nDim+2;

  /*--- Compute Cvvs ---*/
  const auto& cvve_eq = ComputeSpeciesCvVibEle(T);

  /*--- Loop through species ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[nEv][iVar] += rhos[iSpecies]/taus[iSpecies]*(cvve_eq[iSpecies]*dTdU[iVar]-cvve[iSpecies]*dTvedU[iVar]);//TODO*Volume;
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_jacobian[nEv][iSpecies] += (eve_eq[iSpecies]-eve[iSpecies])/taus[iSpecies];//TODO *Volume;
}

vector<su2double>& CSU2TCLib::ComputeSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves){

  vector<su2double> cvtrs;
  //TODO: ADD Electrons?
  cvtrs = GetSpeciesCvTraRot();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    eves[iSpecies] = val_eves[iSpecies];
    hs[iSpecies] = Ru/MolarMass[iSpecies]*val_T + cvtrs[iSpecies]*val_T + Enthalpy_Formation[iSpecies] + eves[iSpecies];
  }

  return hs;

}

vector<su2double>& CSU2TCLib::GetDiffusionCoeff(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
   DiffusionCoeffWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
   DiffusionCoeffGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
   DiffusionCoeffWBE();

  return DiffusionCoeff;

}

su2double CSU2TCLib::GetViscosity(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
    ViscosityWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
    ViscosityGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
    ViscositySuth();

  return Mu;

}

vector<su2double>& CSU2TCLib::GetThermalConductivities(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
    ThermalConductivitiesWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
    ThermalConductivitiesGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
    ThermalConductivitiesSuth();

  if (single_temperature) {
    ThermalConductivities[0] += ThermalConductivities[1];
    ThermalConductivities[1] = 0.0;
  }

  return ThermalConductivities;

}

void CSU2TCLib::DiffusionCoeffWBE(){

  /*--- Calculate species mole fraction ---*/
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;

  /*--- Calculate mixture molar mass (kg/mol) ---*/
  // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
  su2double M = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    M += MolarMass[iSpecies]*MolarFracWBE[iSpecies];
  M = M*1E-3;

  /*---+++                  +++---*/
  /*--- Diffusion coefficients ---*/
  /*---+++                  +++---*/
  /*--- Solve for binary diffusion coefficients ---*/
  // Note: Dij = Dji, so only loop through req'd indices
  // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
  su2activematrix Dij;
  Dij.resize(nSpecies, nSpecies) = su2double(0.0);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    const su2double Mi = MolarMass[iSpecies]*1E-3;
    for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj = MolarMass[jSpecies]*1E-3;

      /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/

      /*--- If collisions between electrons/ion, used Coloumb potentials ---*/
      bool coulomb = false;
      if (abs(Omega11(iSpecies, jSpecies, 0)) == 1.0 && ionization) coulomb = true;

      // Used Tve for electron collisions
      const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T; 

      /*--- Compute the collisional cross section (omega_ij) ---*/
      const su2double Omega_ij = ComputeCollisionCrossSection(iSpecies, jSpecies, T_col, true, coulomb) / PI_NUMBER;

      /*--- Calculate and populate diffusion coefficients ---*/
      Dij(iSpecies,jSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
      Dij(jSpecies,iSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
    }
  }

  /*--- Calculate species-mixture diffusion coefficient --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DiffusionCoeff[iSpecies] = 0.0;
    su2double denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) {
        denom += MolarFracWBE[jSpecies]/Dij(iSpecies,jSpecies);
      }
    }

    if (nSpecies==1) DiffusionCoeff[0] = 0;
    else DiffusionCoeff[iSpecies] = (1-MolarFracWBE[iSpecies])/denom;
  }
}

void CSU2TCLib::ViscosityWBE(){

  /*--- Calculate species mole fraction ---*/
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                             Blottner[iSpecies][1])*log(T) +
                             Blottner[iSpecies][2]);

  /*--- Determine species 'phi' value for Blottner model ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    phis[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(MolarMass[jSpecies]/MolarMass[iSpecies], 0.25);
      const su2double tmp2 = sqrt(8.0*(1.0+MolarMass[iSpecies]/MolarMass[jSpecies]));
      phis[iSpecies] += MolarFracWBE[jSpecies]*tmp1*tmp1/tmp2;
    }
  }

  /*--- Calculate mixture laminar viscosity ---*/
  Mu = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    Mu += MolarFracWBE[iSpecies]*mus[iSpecies]/phis[iSpecies];
  }
}

void CSU2TCLib::ThermalConductivitiesWBE(){

  vector<su2double> ks, kves;

  ks.resize(nSpecies,0.0);
  kves.resize(nSpecies,0.0);

  Cvves = ComputeSpeciesCvVibEle(Tve);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + RotationModes[iSpecies]/2.0)*Ru/MolarMass[iSpecies];
    kves[iSpecies] = mus[iSpecies]*Cvves[iSpecies];
  }

  /*--- Calculate mixture tr & ve conductivities ---*/
  ThermalCond_tr = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ThermalCond_tr += MolarFracWBE[iSpecies]*ks[iSpecies]/phis[iSpecies];
    ThermalCond_ve += MolarFracWBE[iSpecies]*kves[iSpecies]/phis[iSpecies];
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;
}

su2double CSU2TCLib::ComputeCollisionCrossSection(unsigned iSpecies, unsigned jSpecies, su2double T, bool d1, bool coulomb) {

  const su2double pi = PI_NUMBER;
  const su2double Na = AVOGAD_CONSTANT;
  
  if (coulomb) {

    const su2double e_cgs = FUND_ELEC_CHARGE_CGS; // CGS unit of fundamental electric charge 
    const su2double kb_cgs = BOLTZMANN_CONSTANT * 1E7; // CGS unit of Boltzmann Constant 
    const su2double ne_cgs = Na * rhos[0] / MolarMass[0] * 1E-6; // CGS unit of electron number density
        
    const su2double debyeLength = sqrt(kb_cgs * T / 4 / pi / ne_cgs / pow(e_cgs,2));
    const su2double T_star = debyeLength / (pow(e_cgs,2) / (kb_cgs * T));

    /*--- Compute the collisionion cross section ---*/
    // Note: Omega11 is used for diffusion, viscosity, translational, internal, and reaction components of
    //       thermal conductivity
    //       Omega22 is used for viscosity and translational components of thermal conductivity

    if (Omega11(iSpecies, jSpecies, 0) == 1.0 && d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D1_a*T_star*(1 - C1_a * exp(-c1_a * T_star))+1);
    } if (Omega11(iSpecies, jSpecies, 0) == -1.0 && d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D1_r*T_star*(1 - C1_r * exp(-c1_r * T_star))+1);
    } else if (Omega22(iSpecies, jSpecies, 0) == 1.0 && !d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D2_a*T_star*(1 - C2_a * exp(-c2_a * T_star))+1);
    } else {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D2_r*T_star*(1 - C2_r * exp(-c2_r * T_star))+1);
    }

  } else {
    if (d1) {
      return 1E-20 * Omega11(iSpecies,jSpecies,3) * pow(T, Omega11(iSpecies,jSpecies,0)*log(T)*log(T) + Omega11(iSpecies,jSpecies,1)*log(T) + Omega11(iSpecies,jSpecies,2));
    }       return 1E-20 * Omega22(iSpecies,jSpecies,3) * pow(T, Omega22(iSpecies,jSpecies,0)*log(T)*log(T) + Omega22(iSpecies,jSpecies,1)*log(T) + Omega22(iSpecies,jSpecies,2));
   
  }
}

su2double CSU2TCLib::ComputeCollisionDelta(unsigned iSpecies, unsigned jSpecies, su2double Mi, su2double Mj, su2double T, bool d1) {

  bool coulomb = false;
  if (abs(Omega11(iSpecies, jSpecies, 0)) == 1.0 && ionization) {
    coulomb = true;
  } 

  const su2double Omega_ij = ComputeCollisionCrossSection(iSpecies, jSpecies, T, d1, coulomb);
  const su2double pi = PI_NUMBER;
  su2double delta = 0.0;

  if (d1) {
    delta = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij; // d1_ij
  } else {
    delta = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij; // d2_ij
  }
  return fmin(delta, 1E16);
}

void CSU2TCLib::DiffusionCoeffGY(){

  /*--- Calculate mixture gas constant ---*/
  su2double gam_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gam_t += rhos[iSpecies] / (Density*MolarMass[iSpecies]);
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);
    su2double denom = 0.0;

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) { 

        const su2double Mj    = (MolarMass[jSpecies] + EPS);
        const su2double gam_j = rhos[jSpecies] / (Density*Mj);

        const su2double kb = BOLTZMANN_CONSTANT;

        const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T; 

        su2double d1_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, true);

        const su2double D_ij = kb*T_col/(Pressure*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    /*--- Calculate species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = (gam_t*gam_t*Mi*(1-Mi*gam_i) / denom);
  }
}

void CSU2TCLib::ViscosityGY(){

  const su2double Na = AVOGAD_CONSTANT;
  Mu = 0.0;

  /*--- Mixture viscosity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    su2double denom = 0.0;

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj    = (MolarMass[jSpecies] + EPS);
      const su2double gam_j = rhos[jSpecies] / (Density*Mj);

      const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T; 

      su2double d2_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, false);

      denom += gam_j*d2_ij;
    }
    /*--- Calculate species laminar viscosity ---*/
    Mu += (Mi/Na * gam_i) / denom;
  }
}

void CSU2TCLib::ThermalConductivitiesGY(){

  const su2double Na   = AVOGAD_CONSTANT;
  const su2double kb   = BOLTZMANN_CONSTANT;

  /*--- Mixture vibrational-electronic specific heat ---*/
  const auto Cvves = ComputeSpeciesCvVibEle(Tve);
  su2double rhoCvve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoCvve += rhos[iSpecies]*Cvves[iSpecies];
  const su2double Cvve = rhoCvve/Density;

  /*--- Calculate mixture gas constant ---*/
  su2double R = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    R += Ru / MolarMass[iSpecies] * rhos[iSpecies]/Density;
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  su2double ThermalCond_tr = 0.0;
  su2double ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double mi    = Mi/Na;
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);
    su2double denom_t = 0.0;
    su2double denom_r = 0.0;
    su2double denom_re = 0.0;

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj    = (MolarMass[jSpecies] + EPS);
      const su2double mj    = Mj/Na;
      const su2double gam_j = rhos[iSpecies] / (Density*Mj);
      const su2double a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));

      const su2double T_col = ((iSpecies == 0 && ionization) || (jSpecies == 0 && ionization)) ? Tve : T; 

      su2double d1_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, true);
      su2double d2_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, false);

      if (jSpecies == 0 && ionization) { denom_t += 3.54*gam_j*d2_ij; }
      else { denom_t += a_ij*gam_j*d2_ij; }

      denom_r += gam_j*d1_ij;
      denom_re += gam_j*d2_ij;
    }

    /*--- Prevent divide by 0 ---*/
    if (denom_t <= 0.0) denom_t = EPS;
    if (denom_r <= 0.0) denom_r = EPS;
    if (denom_re <= 0.0) denom_re = EPS;

    /*--- Translational contribution to thermal conductivity ---*/
    if (!ionization || iSpecies != 0) ThermalCond_tr += ((15.0/4.0)*kb*gam_i/denom_t);

    /*--- Rotational contribution to thermal conductivity ---*/
    if (RotationModes[iSpecies] != 0.0) ThermalCond_tr += (kb*gam_i/denom_r);

    /*--- Vibrational-electronic contribution to thermal conductivity ---*/
    if ((!ionization || iSpecies != 0) && RotationModes[iSpecies] != 0.0) ThermalCond_ve += (kb*Cvve/R*gam_i / denom_r);
    
    if (ionization && iSpecies == 0) ThermalCond_ve += ((15.0/4.0)*kb*gam_i/(1.45*denom_re));
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;
}

void CSU2TCLib::ViscositySuth(){

  su2double T_nd = T / T_ref_suth;

  /*--- Calculate mixture laminar viscosity ---*/
  Mu = mu_ref[0] * T_nd * sqrt(T_nd) * ((T_ref_suth + Sm_ref[0]) / (T + Sm_ref[0]));
}

void CSU2TCLib::ThermalConductivitiesSuth(){

  /*--- Compute mixture quantities ---*/
  su2double mass = 0.0, rho = 0.0;
  for (unsigned short ii=0; ii<nSpecies; ii++) rho  += rhos[ii];
  for (unsigned short ii=0; ii<nSpecies; ii++) mass += rhos[ii]/rho*MolarMass[ii];

  su2double Cvtr = ComputerhoCvtr()/rho;
  su2double Cvve = ComputerhoCvve()/rho;

  /*--- Compute simple Kve scaling factor ---*/
  su2double scl  = Cvve/Cvtr;

  /*--- Compute k's using Sutherland's law ---*/
  su2double T_nd = T / T_ref_suth;
  su2double k = k_ref[0] * T_nd * sqrt(T_nd) * ((T_ref_suth + Sk_ref[0]) / (T + Sk_ref[0]));
  su2double kve = scl*k;

  ThermalConductivities[0] = k;
  ThermalConductivities[1] = kve;
}

vector<su2double>& CSU2TCLib::ComputeTemperatures(vector<su2double>& val_rhos, su2double rhoE, su2double rhoEve, su2double rhoEvel, su2double Tve_old) {

  rhos = val_rhos;

  /*----------Translational temperature----------*/
  su2double rhoE_f   = 0.0;
  su2double rhoE_ref = 0.0;
  su2double rhoCvtr  = 0.0;
  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
    rhoCvtr  += rhos[iSpecies] * Cvtrs[iSpecies];
    rhoE_ref += rhos[iSpecies] * Cvtrs[iSpecies] * Ref_Temperature[iSpecies];
    rhoE_f   += rhos[iSpecies] * (Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies]);
  }

  T = (rhoE - rhoEve - rhoE_f + rhoE_ref - rhoEvel) / rhoCvtr;

  /*--- Set temperature clipping values ---*/
  const su2double Tmin   = 50.0; const su2double Tmax   = 8E4;
  const su2double Tvemin = 50.0; const su2double Tvemax = 8E4;
  su2double Tve_o  = 50.0; su2double Tve2  = 8E4;

  /* Determine if the temperature lies within the acceptable range */
  if (Tve_old < 1) Tve_old = T;                           //For first fluid iteration
  if (T < Tmin) T = Tmin;  else if (T > Tmax) T = Tmax;
  if (Tve_old<Tvemin) Tve_old = Tvemin; else if (Tve_old>Tvemax) Tve_old = Tvemax;

  /*--- Set vibrational temperature algorithm parameters ---*/
  const su2double NRtol         = 1.0E-6;    // Tolerance for the Newton-Raphson method
  const su2double Btol          = 1.0E-6;    // Tolerance for the Bisection method
  const unsigned short maxBIter = 100;        // Maximum Bisection method iterations
  const unsigned short maxNIter = 100;        // Maximum Newton-Raphson iterations
  const su2double scale         = 0.9;       // Scaling factor for Newton-Raphson step

  /*--- Execute a Newton-Raphson root-finding method for Tve ---*/
  //Initialize solution
  Tve = Tve_old;

  bool Bconvg = false;
  bool NRconvg = false;
  su2double rhoEve_t = 0.0, rhoCvve = 0.0;

  /*--- Newton-Raphson Method --*/
  for (unsigned short iIter = 0; iIter < maxNIter; iIter++) {
    rhoEve_t = rhoCvve = 0.0;
    const auto& val_eves  = ComputeSpeciesEve(Tve);
    const auto& val_cvves = ComputeSpeciesCvVibEle(Tve);

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
      rhoEve_t += rhos[iSpecies] * val_eves[iSpecies];
      rhoCvve += rhos[iSpecies] * val_cvves[iSpecies];
    }

    /*--- Find the roots ---*/
    su2double f  = rhoEve - rhoEve_t;
    su2double df = -rhoCvve;
    Tve2 = Tve - (f/df)*scale;

    /*--- Check for convergence ---*/
    if ((fabs(Tve2-Tve) < NRtol) && (Tve > Tvemin) && (Tve < Tvemax)) {
      NRconvg = true;
      Tve = Tve2;
      break;
    }       Tve = Tve2;
   
  }

  // If the Newton-Raphson method has converged, assign the value of Tve.
  // Otherwise, execute a bisection root-finding method
  Tve_o = Tvemin; Tve2 = Tvemax;
  if (!NRconvg) {
    for (unsigned short iIter = 0; iIter < maxBIter; iIter++) {
      Tve      = (Tve_o+Tve2)/2.0;
      const auto& val_eves = ComputeSpeciesEve(Tve);
      rhoEve_t = 0.0;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) rhoEve_t += rhos[iSpecies] * val_eves[iSpecies];
      if (fabs(rhoEve_t - rhoEve) < Btol) {
        Bconvg = true;
        break;
      }         if (rhoEve_t > rhoEve) Tve2 = Tve;
        else                  Tve_o = Tve;
     
    }
  }

  // If absolutely no convergence, then assign to the TR temperature
  if (!NRconvg && !Bconvg ) {
    Tve = T;
  }

  if (single_temperature) {
    Tve = T;
  }

  temperatures[0] = T;
  temperatures[1] = Tve;

  return temperatures;
}

void CSU2TCLib::GetChemistryEquilConstants(unsigned short iReaction){

  if (custom_chemistry_tables) {
    if (iReaction >= CustomKeqCoeff.size()) {
      SU2_MPI::Error("CONFIG ERROR: Missing custom Keq coefficients for requested reaction index.", CURRENT_FUNCTION);
    }
    for (unsigned short iRow = 0; iRow < 6; ++iRow) {
      for (unsigned short iCol = 0; iCol < 5; ++iCol) {
        RxnConstantTable(iRow, iCol) = CustomKeqCoeff[iReaction][iCol];
      }
    }
    return;
  }

  if (gas_model == "O2"){
    // THESE ARE UNUSED.  SHOULD WE KEEP????  Good for future?
    //O2 + M -> 2O + M
    RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
    RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
    RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
    RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
    RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
    RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

  } else if (gas_model == "N2"){

    //N2 + M -> 2N + M
    RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
    RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
    RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
    RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
    RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
    RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

  } else if (gas_model == "AIR-5"){

    if (iReaction <= 4) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 4 && iReaction <= 9) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 9 && iReaction <= 14) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 15) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 16) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;
    }

  } else if (gas_model == "AIR-7"){

    if (iReaction <= 5) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 5 && iReaction <= 11) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 11 && iReaction <= 17) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 18) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 19) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;

    } else if (iReaction == 20) {

      //N + O -> NO+ + e-
      RxnConstantTable(0,0) = -2.1852;   RxnConstantTable(0,1) = -6.6709; RxnConstantTable(0,2) = -4.2968; RxnConstantTable(0,3) = -2.2175; RxnConstantTable(0,4) = -0.050748;
      RxnConstantTable(1,0) = -1.0276;   RxnConstantTable(1,1) = -7.1278; RxnConstantTable(1,2) = -2.637;  RxnConstantTable(1,3) = -2.95;   RxnConstantTable(1,4) = -0.0021;
      RxnConstantTable(2,0) = -0.65871;  RxnConstantTable(2,1) = -7.2742; RxnConstantTable(2,2) = -2.1096; RxnConstantTable(2,3) = -3.1823; RxnConstantTable(2,4) = 0.01331;
      RxnConstantTable(3,0) = -0.57924;  RxnConstantTable(3,1) = -7.3079; RxnConstantTable(3,2) = -1.9999; RxnConstantTable(3,3) = -3.2294; RxnConstantTable(3,4) = 0.016382;
      RxnConstantTable(4,0) = -0.53538;  RxnConstantTable(4,1) = -7.3252; RxnConstantTable(4,2) = -1.937;  RxnConstantTable(4,3) = -3.2572; RxnConstantTable(4,4) = 0.01823;
      RxnConstantTable(5,0) = -0.52801;  RxnConstantTable(5,1) = -7.3281; RxnConstantTable(5,2) = -1.9264; RxnConstantTable(5,3) = -3.2618; RxnConstantTable(5,4) = 0.01854;

    } else if (iReaction == 21) {

      //N2 + e -> N + N + e
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978; RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.3897;  RxnConstantTable(1,2) = 2.0617; RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.5732;  RxnConstantTable(2,2) = 1.3923; RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543;
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.6061;  RxnConstantTable(3,2) = 1.2993; RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.6291;  RxnConstantTable(4,2) = 1.2153; RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.009444;
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.6291;  RxnConstantTable(5,2) = 1.2153; RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.009444;
    }
  }
}
