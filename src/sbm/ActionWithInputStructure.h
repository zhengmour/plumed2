/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2022 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_sbm_ActionWithInputStructure_H
#define __PLUMED_sbm_ActionWithInputStructure_H

#include "adjmat/AdjacencyMatrixBase.h"
#include "adjmat/AdjacencyMatrixVessel.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/GenericMolInfo.h"
#include "tools/PDB.h"

#include <vector>
#include <string>

namespace PLMD{
namespace sbm{

class PDB;
class GenericMolInfo;
class AdjacencyMatrixVessel;

class Bond {
  int v;
  int w;
  double distance;
public:
  Bond(int v, int w, double distance) :
    v(v), w(w), distance(distance) { }
  double getDistance() { return distance; }
  int either() { return v; }
  int other(int vertex) { 
    if      (vertex == v) return w;
    else if (vertex == w) return v;
    else throw "Inconsistent bond";
  }
  int compareTo(Bond* that) {
    if      (this->distance < that->distance) return -1;
    else if (this->distance > that->distance) return +1;
    else                                      return  0;
  }
};

class ActionWithInputStructure : public ActionAtomistic
{
protected:
  ForwardDecl<AdjacencyMatrixVessel> mymatrix_pwd;
/// The vessel that holds the adjacency matrix
  AdjacencyMatrixVessel& mymatrix=*mymatrix_pwd;

  ForwardDecl<PDB> pdb_pwd;
/// a pdb file containing the reference structure
  PDB& pdb=*pdb_pwd;
/// The type of molecule in the pdb
  std::string mytype;
/// The name of the reference pdb file
  std::string reference;
/// The backbone that was read in from the pdb file
  std::vector< std::vector<AtomNumber> > read_backbone;
/// Structure in pdb file is whole 
  bool iswhole_;

/// The vector that contains contact relations
  std::vector<Bond> nativeBonds;
  std::vector<Bond> nonativeBonds;
/// The coefficient of each contact relation
  std::vector<double> nativeEplisons;   // native contact
  std::vector<double> nonativeEplisons; // non-native contact
  std::vector<double> excvolEplisons;   // excluded volume
  std::vector<double> excvolSigmas;
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionWithInputStructure(const ActionOptions&);
  ~ActionWithInputStructure();

/// Retrieve the vessel that holds the adjacency matrix
  AdjacencyMatrixVessel* getAdjacencyVessel() const;
/// Retrieve the vector that contains native contact informations
  void getNativeBonds(std::vector<Bond>& bonds) const;
/// Retrieve the vector that contains 
  void getNonativeBonds(std::vector<Bond>& bonds) const;

/// read in some atoms
  bool parseStructureAtomList(const std::string& key, const int& num, std::vector<AtomNumber>& t);
/// read contact infomations
  bool readContactInfo(const std::string& key, std::vector<Bond>& bonds, std::vector<double>& eplisons);
/// read excluded volume 
  bool readExcludedVolumeInfo(const std::string& key, std::vector<double>& eplisons, std::vector<double>& sigmas);

/// Retrieve the information containing in pdb file
  void getBackbone( std::vector<std::string>& resstrings, const std::string& fortype, std::vector< std::vector<AtomNumber> >& backbone );
  std::string getAtomName(AtomNumber a)const;
  bool checkForAtom(AtomNumber a)const;
  Vector getPosition(AtomNumber a) const;
  bool isWhole();
  bool isPeriodic() { return false; }
/// No loop over tasks for ActionWithInputStructure
  double compute() { plumed_error(); }
};

}
}


#endif ! __PLUMED_sbm_ActionWithInputMatrix_H
