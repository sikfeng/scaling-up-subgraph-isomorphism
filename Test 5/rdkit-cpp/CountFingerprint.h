#ifndef _RD_COUNTFINGERPRINTS_H_
#define _RD_COUNTFINGERPRINTS_H_

#include <vector>
#include <cstdint>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>

//class RDKit::ROMol;

//! \brief Generates a topological (Daylight like) fingerprint for a molecule
//!        using an alternate (faster) hashing algorithm
/*!

  \param mol:          the molecule to be fingerprinted
  \param minPath:      the minimum path length (in bonds) to be included
  \param maxPath:      the minimum path length (in bonds) to be included
  \param fpSize:       the size of the fingerprint
  \param nBitsPerHash: the number of bits to be set by each path
  \param useHs:        toggles inclusion of Hs in paths (if the molecule has
  explicit Hs)
  \param tgtDensity:   if the generated fingerprint is below this density, it
  will
                       be folded until the density is reached.
  \param minSize:      the minimum size to which the fingerprint will be
                       folded
  \param branchedPaths: toggles generation of branched subgraphs, not just
  linear paths
  \param useBondOrders: toggles inclusion of bond orders in the path hashes
  \param atomInvariants: a vector of atom invariants to use while hashing the
  paths
  \param fromAtoms:    only paths starting at these atoms will be included
  \param atomBits:     used to return the bits that each atom is involved in
                       (should be at least \c mol.numAtoms long)

  \return the molecular fingerprint, as an ExplicitBitVect

  <b>Notes:</b>
    - the caller is responsible for <tt>delete</tt>ing the result

*/

std::vector<uint32_t> CountFingerprintMol(
    const RDKit::ROMol &mol, unsigned int minPath = 1, unsigned int maxPath = 7,
    unsigned int fpSize = 2048, unsigned int nBitsPerHash = 2,
    bool useHs = true, double tgtDensity = 0.0, unsigned int minSize = 128,
    bool branchedPaths = true, bool useBondOrder = true,
    std::vector<std::uint32_t> *atomInvariants = 0,
    const std::vector<std::uint32_t> *fromAtoms = 0);

#endif
