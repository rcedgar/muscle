/////////////////////////////////////////////////////////////////
// EvolutionaryTree.h
//
// Utilities for reading/writing multiple sequence data.
/////////////////////////////////////////////////////////////////

#ifndef EVOLUTIONARYTREE_H
#define EVOLUTIONARYTREE_H

#include <string>
#include <list>
#include <stdio.h>
#include "multisequence.h"
#include "sequence.h"

using namespace std;

/////////////////////////////////////////////////////////////////
// TreeNode
//
// The fundamental unit for representing an alignment tree.  The
// guide tree is represented as a binary tree.
/////////////////////////////////////////////////////////////////

class TreeNode {
  int sequenceIndex;                  // sequence label
  TreeNode *left, *right, *parent;    // pointers to left, right children

  /////////////////////////////////////////////////////////////////
  // TreeNode::PrintNode()
  //
  // Internal routine used to print out the sequence comments
  // associated with the evolutionary tree, using a hierarchical
  // parenthesized format.
  /////////////////////////////////////////////////////////////////

  void PrintNode (ostream &outfile, const MultiSequence *sequences) const;
 public:

  /////////////////////////////////////////////////////////////////
  // TreeNode::TreeNode()
  //
  // Constructor for a tree node.  Note that sequenceIndex = -1
  // implies that the current node is not a leaf in the tree.
  /////////////////////////////////////////////////////////////////

  TreeNode (int sequenceIndex) : sequenceIndex (sequenceIndex),
    left (NULL), right (NULL), parent (NULL) {
    assert (sequenceIndex >= -1);
  }

  /////////////////////////////////////////////////////////////////
  // TreeNode::~TreeNode()
  //
  // Destructor for a tree node.  Recursively deletes all children.
  /////////////////////////////////////////////////////////////////

  ~TreeNode (){
    if (left){ delete left; left = NULL; }
    if (right){ delete right; right = NULL; }
    parent = NULL;
  }


  // getters
  int GetSequenceLabel () const { return sequenceIndex; }
  TreeNode *GetLeftChild () const { return left; }
  TreeNode *GetRightChild () const { return right; }
  TreeNode *GetParent () const { return parent; }

  // setters
  void SetSequenceLabel (int sequenceIndex){ this->sequenceIndex = sequenceIndex; assert (sequenceIndex >= -1); }
  void SetLeftChild (TreeNode *left){ this->left = left; }
  void SetRightChild (TreeNode *right){ this->right = right; }
  void SetParent (TreeNode *parent){ this->parent = parent; }

  /////////////////////////////////////////////////////////////////
  // TreeNode::ComputeTree()
  //
  // Routine used to compute an evolutionary tree based on the
  // given distance matrix.  We assume the distance matrix has the
  // form, distMatrix[i][j] = expected accuracy of aligning i with j.
  /////////////////////////////////////////////////////////////////

  static TreeNode *ComputeTree (const vector<vector<float> > &distMatrix);
  /////////////////////////////////////////////////////////////////
  // TreeNode::Print()
  //
  // Print out the subtree associated with this node in a
  // parenthesized representation.
  /////////////////////////////////////////////////////////////////

  void Print (ostream &outfile, const MultiSequence *sequences) const;
};

#endif
