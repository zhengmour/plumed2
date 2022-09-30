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
#include "ActionWithInputStructure.h"
#include "vesselbase/ActionWithVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace sbm{

void ActionWithInputStructure::registerKeywords(Keywords& keys) {
  ActionSetup::registerKeywords(keys);
  keys.add("compulsory", "STRUCTURE", "a file in pdb format containing a reference structure. "
           "This is used to defines the atoms in the various residues, chains, etc . "
           "For more details on the PDB file format visit http://www.wwpdb.org/docs.html");
  keys.add("compulsory", "MOLTYPE", "zeolite", "what kind of molecule is contained in pdb file.");
  keys.add("compulsory", "CONTACT", "which contains native contact and excluded volume, "
           "additional excluded volume has a large effect on the entropy of the structure "
           "where most contacts are formed.");
  keys.add("optional", "NON_NATIVE", "which contains non-native contacts.");
  keys.addFlag("WHOLE", false, "The reference structure is whole, i.e. not broken by PBC");
}

ActionWithInputStructure::~ActionWithInputStructure() { }

ActionWithInputStructure::ActionWithInputStructure(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  iswhole_(false)
{
  // Read what is containd in the pdb file
  parse("MOLTYPE", mytype);

  // check if whole
  parseFlag("WHOLE", iswhole_);

  auto* moldat=plumed.getActionSet().selectLatest<ActionWithInputStructure*>(this);
  if( moldat ) log<<"  overriding last MOLINFO with label " << moldat->getLabel()<<"\n";

  // Read filename of pdb file
  parse("STRUCTURE", reference);

  if( ! pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()))plumed_merror("missing input file " + reference );

  std::vector<std::string> chains; pdb.getChainNames( chains );
  log.printf("  pdb file named %s contains %u chains \n",reference.c_str(), static_cast<unsigned>(chains.size()));
  for(unsigned i=0; i<chains.size(); ++i) {
    unsigned start,end; std::string errmsg;
    pdb.getResidueRange( chains[i], start, end, errmsg );
    if( errmsg.length()!=0 ) error( errmsg );
    AtomNumber astart,aend;
    pdb.getAtomRange( chains[i], astart, aend, errmsg );
    if( errmsg.length()!=0 ) error( errmsg );
    log.printf("  chain named %s contains residues %u to %u and atoms %u to %u \n",chains[i].c_str(),start,end,astart.serial(),aend.serial());
  }

  // Read coefficients of contact relations
  parse("CONTACT", contactstr);

  // Read coefficients of non-native relations
  parse("NON_NATIVE", nonativestr);  



}

}
}


