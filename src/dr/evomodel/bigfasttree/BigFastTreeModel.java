/*
 * BigFastTreeModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.bigfasttree;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Iterator;

/**
 * A model component for trees. This version uses flat arrays for nodes and heights to maximize
 * speed (particularly storing/restoring) on very large trees.
 *
 * @author Andrew Rambaut
 */
public class BigFastTreeModel extends TreeModel {

    //
    // Public stuff
    //

    public static final String BIG_FAST_TREE_MODEL = "bigFastTreeModel";

    public BigFastTreeModel(Tree tree) {
        this(BIG_FAST_TREE_MODEL, tree, false, false);
    }

    public BigFastTreeModel(String name, Tree tree) { this(name, tree, false, false); }

    /* New constructor that copies the attributes of Tree tree into the new TreeModel
      * Useful for constructing a TreeModel from a NEXUS file entry
      */

    public BigFastTreeModel(String name, Tree tree, boolean fixHeights, boolean fixTree) {

        super(name, !fixTree);
        setId(name);

        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree);
        binaryTree.resolveTree();

        // adjust the heights to be compatible with the tip dates and perturb
        // any zero branches.
        if (!fixHeights) {
            MutableTree.Utils.correctHeightsForTips(binaryTree);
        }

        internalNodeCount = binaryTree.getInternalNodeCount();
        externalNodeCount = binaryTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new NodeRef[nodeCount];

        edges = new int[nodeCount * 3];
        storedEdges = new int[nodeCount * 3];

        heights = new double[nodeCount];
        storedHeights = new double[nodeCount];

        boolean done = false;
        NodeRef node = binaryTree.getRoot();
        do {
            node = TreeUtils.postorderSuccessor(binaryTree, node);

            int number = node.getNumber();

            if (binaryTree.isExternal(node)) {
                nodes[number] = new Node(number, binaryTree.getNodeTaxon(node));
                edges[(number * 3)] = binaryTree.getParent(node).getNumber(); // parent
                edges[(number * 3) + 1] = -1; // child 1
                edges[(number * 3) + 2] = -1; // child 2
            } else {
                if (binaryTree.isRoot(node)) {
                    root = number;
                    done = true;
                    edges[(number * 3)] = -1;
                } else {
                    edges[(number * 3)] = binaryTree.getParent(node).getNumber(); // parent
                }
                nodes[number] = new Node(number);

                edges[(number * 3) + 1] = binaryTree.getChild(node, 0).getNumber(); // child 1
                edges[(number * 3) + 2] = binaryTree.getChild(node, 1).getNumber(); // child 2;
            }
            heights[number] = binaryTree.getNodeHeight(node);

        } while (!done);
    }

    protected void copyTopology(Tree tree){
        // get a rooted version of the tree to clone
        FlexibleTree binaryTree = new FlexibleTree(tree);
        binaryTree.resolveTree();

        boolean done = false;
        NodeRef node = binaryTree.getRoot();
        do {
            node = TreeUtils.postorderSuccessor(binaryTree, node);

            int number = node.getNumber();

            if (binaryTree.isExternal(node)) {
                nodes[number] = new Node(number, binaryTree.getNodeTaxon(node));
                edges[(number * 3)] = binaryTree.getParent(node).getNumber(); // parent
                edges[(number * 3) + 1] = -1; // child 1
                edges[(number * 3) + 2] = -1; // child 2
            } else {
                if (binaryTree.isRoot(node)) {
                    root = number;
                    done = true;
                    edges[(number * 3)] = -1;
                } else {
                    edges[(number * 3)] = binaryTree.getParent(node).getNumber(); // parent
                }
                nodes[number] = new Node(number);

                edges[(number * 3) + 1] = binaryTree.getChild(node, 0).getNumber(); // child 1
                edges[(number * 3) + 2] = binaryTree.getChild(node, 1).getNumber(); // child 2;
            }
            heights[number] = binaryTree.getNodeHeight(node);

        } while (!done);
    };


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no submodels so nothing to do
    }

    /**
     * Called when a parameter changes.
     */
    @Override
    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no parameters so nothing to do
    }

    @Override
    public boolean inTreeEdit() {
        return inEdit;
    }

    // *****************************************************************
    // Interface Tree
    // *****************************************************************

    /**
     * @return a count of the number of nodes (internal + external) in this
     *         tree.
     */
    @Override
    public int getNodeCount() {
        return nodeCount;
    }

    @Override
    public double getNodeHeight(NodeRef node) {
        return heights[node.getNumber()];
    }

    public double getNodeHeightUpper(NodeRef node) {
        if (isRoot(node)) {
            return Double.POSITIVE_INFINITY;
        }

        return getNodeHeight(getParent(node));
    }

    public double getNodeHeightLower(NodeRef node) {
        if (isExternal(node)) {
            return 0.0;
        }
        return Math.max(
                getNodeHeight(getChild(node, 0)),
                getNodeHeight(getChild(node, 1)));
    }


    /**
     * @param node
     * @return the rate parameter associated with this node.
     */
    @Override
    public double getNodeRate(NodeRef node) {
//        //mimics default tree that does not have node rates
            return 1.0;
//        throw new UnsupportedOperationException("getNodeRate not available in BigFastTreeModel");
    }

    @Override
    public Object getNodeAttribute(NodeRef node, String name) {
        throw new UnsupportedOperationException("getNodeAttribute not available in BigFastTreeModel");
    }

    @Override
    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new UnsupportedOperationException("getNodeAttributeNames not available in BigFastTreeModel");
    }

    @Override
    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new UnsupportedOperationException("getMultivariateNodeTrait not available in BigFastTreeModel");
    }

    @Override
    public Taxon getNodeTaxon(NodeRef node) {
        return ((Node) node).taxon;
    }

    @Override
    public boolean isExternal(NodeRef node) {
        return node.getNumber() < externalNodeCount;
    }

    @Override
    public boolean isRoot(NodeRef node) {
        return node.getNumber() == root;
    }

    @Override
    public int getChildCount(NodeRef node) {
        if(isExternal(node)){
            return 0;
        }
        int kids=0;
        for (int i = 0; i < 2; i++) {
            if(edges[(node.getNumber() * 3) + i + 1]>-1){
                kids+=1;
            };
        }

        return kids;

    }

    @Override
    public NodeRef getChild(NodeRef node, int i) {
        return nodes[getChild(node.getNumber(), i)];
    }

    @Override
    public NodeRef getParent(NodeRef node) {
        if (isRoot(node)) {
            return null;
        }
        int parentIndex = getParent(node.getNumber());
        return parentIndex == -1 ? null : nodes[parentIndex];
    }

    @Override
    public NodeRef getExternalNode(int i) {
        return nodes[i];
    }

    @Override
    public NodeRef getInternalNode(int i) {
        return nodes[i + externalNodeCount];
    }

    @Override
    public NodeRef getNode(int i) {
        return nodes[i];
    }

    @Override
    public NodeRef[] getNodes() {
        return nodes;
    }

    /**
     * Returns the number of external nodes.
     */
    @Override
    public int getExternalNodeCount() {
        return externalNodeCount;
    }

    /**
     * Returns the ith internal node.
     */
    @Override
    public int getInternalNodeCount() {
        return internalNodeCount;
    }

    /**
     * Returns the root node of this tree.
     */
    @Override
    public NodeRef getRoot() {
        return nodes[root];
    }

    @Override
    public boolean isTipDateSampled() {
        return false;
    }

    private int getParent(int nodeNumber) {
        return edges[(nodeNumber * 3)];
    }

    private int getChild(int nodeNumber, int i) {
        assert i == 0 || i == 1;
        return edges[(nodeNumber * 3) + i + 1];
    }

    private void setParent(int nodeNumber, int parentNumber) {
         edges[(nodeNumber * 3)] = parentNumber;
    }

    private void setChild(int nodeNumber, int i, int childNumber) {
        assert i == 0 || i == 1;
        edges[(nodeNumber * 3) + i + 1] = childNumber;
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    @Override
    public void setRoot(NodeRef newRoot) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        root = newRoot.getNumber();
    }

    @Override
    public void addChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        int parent = p.getNumber();
        int child = c.getNumber();

        if (getChild(parent, 0) == -1) {
            setChild(parent, 0, child);
        } else if (getChild(parent, 1) == -1) {
            setChild(parent, 1, child);
        } else {
            throw new IllegalArgumentException("Node already has two children");
        }
        setParent(child, parent);

        pushTreeChangedEvent(TreeChangedEvent.create(p, false));
        pushTreeChangedEvent(TreeChangedEvent.create(c, false));
    }

    @Override
    public void removeChild(NodeRef p, NodeRef c) {

        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        int parent = p.getNumber();
        int child = c.getNumber();

        if (getChild(parent, 0) == child) {
            setChild(parent, 0, -1);
            //move other child up
            if(getChildCount(p)==1){
                setChild(parent,0,getChild(parent,1));
                setChild(parent, 1, -1);
            }

        } else if (getChild(parent, 1) == child) {
            setChild(parent, 1, -1);
        } else {
            throw new IllegalArgumentException("Child not in node");
        }
        setParent(child, -1);

        pushTreeChangedEvent(TreeChangedEvent.create(p, false));
        pushTreeChangedEvent(TreeChangedEvent.create(c, false));
    }

    @Override
    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new RuntimeException("Unimplemented");
    }


    @Override
    public boolean isTreeValid() {
        for (NodeRef node : nodes) {
            double height = getNodeHeight(node);
            if (height > getNodeHeightUpper(node) || height < getNodeHeightLower(node)) {
                return false;
            }
        }
        return true;
    }

    @Override
    public void setNodeHeight(NodeRef node, double height) {
        heights[node.getNumber()] = height;
        pushTreeChangedEvent(TreeChangedEvent.create(node, true));
    }

    @Override
    public void setNodeHeightQuietly(NodeRef n, double height) {
        heights[n.getNumber()] = height;
    }

    @Override
    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have node rates set");
    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have traits set");
    }

    @Override
    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have traits set");
    }

    @Override
    public void setBranchLength(NodeRef node, double length) {
        throw new UnsupportedOperationException("BigFastTreeModel cannot have branch lengths set");
    }

    @Override
    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new UnsupportedOperationException("BigFastTreeModel does not use NodeAttributes");
    }

    /**
     * Modifies the current tree by adopting the provided collection of edges
     * @param edges Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder Array that contains whether a child node is left or right child
     */
    public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {

        int[] nodeMap = createNodeMap(taxaNames);

        if (this.nodeCount != edges.length) {
            throw new RuntimeException("Incorrect number of edges provided: " + edges.length + " versus " + this.nodeCount + " nodes.");
        }

        //first remove all the child nodes of the internal nodes
        for (int i = this.externalNodeCount; i < this.nodeCount; i++) {
            int childCount = this.getChildCount(nodes[i]);
            for (int j = 0; j < childCount; j++) {
                this.removeChild(nodes[i],this.getChild(nodes[i],0)); // nodes move into the first spot when it's filled
            }
        }

        //start with setting the external node heights
        for (int i = 0; i < this.getExternalNodeCount(); i++) {
            this.setNodeHeight(this.getExternalNode(nodeMap[i]), nodeHeights[i]);
        }
        //set the internal node heights
        for (int i = 0; i < (this.getExternalNodeCount() - 1); i++) {
            //No just restart counting, will fix later on in the code by adding additionalTaxa variable
            this.setNodeHeight(this.getInternalNode(i), nodeHeights[this.getExternalNodeCount() + i]);
        }

        int newRootIndex = -1;
        //now add the parent-child links again to ALL the nodes
        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                //make distinction between external nodes and internal nodes
                if (i < this.getExternalNodeCount()) {
                    //external node
                    this.addChild(this.getNode(edges[i]), this.getExternalNode(nodeMap[i]));
                    System.out.println("external: " + edges[i] + " > " + nodeMap[i]);
                } else {
                    //internal node
                    this.addChild(this.getNode(edges[i]), this.getNode(i));
                    System.out.println("internal: " + edges[i] + " > " + i);
                }
            } else {
                newRootIndex = i;
            }
        }

        //not possible to determine correct ordering of child nodes in the loop where they're being assigned
        //hence perform possible swaps in a separate loop

        for (int i = 0; i < edges.length; i++) {
            if (edges[i] != -1) {
                if(i < this.externalNodeCount) {
                    if (childOrder[i] == 0 && getChild(nodes[edges[i]],0) != nodes[nodeMap[i]]) {
                        //swap child nodes
                        NodeRef node = getNode(edges[i]);

                        NodeRef childOne = getChild(node,0);
                        NodeRef childTwo = getChild(node,1);

                        removeChild(node,childOne);
                        removeChild(node,childTwo);

                        addChild(node,childTwo);
                        addChild(node,childOne);
                    }
                }else{
                    if (childOrder[i] == 0 && getChild(nodes[edges[i]],0) != nodes[i]) {
                        //swap child nodes
                        NodeRef node = getNode(edges[i]);

                        NodeRef childOne = getChild(node,0);
                        NodeRef childTwo = getChild(node,1);

                        removeChild(node,childOne);
                        removeChild(node,childTwo);

                        addChild(node,childTwo);
                        addChild(node,childOne);
                    }
                }
            }

        }

        this.setRoot(nodes[newRootIndex]);
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    /**
     * Store current state
     */
    @Override
    protected void storeState() {
        System.arraycopy(edges, 0, storedEdges, 0, edges.length);
        System.arraycopy(heights, 0, storedHeights, 0, heights.length);

        storedRoot = root;

    }

    /**
     * Restore the stored state
     */
    @Override
    protected void restoreState() {

        int[] tmp = storedEdges;
        storedEdges = edges;
        edges = tmp;

        double[] tmp2 = storedHeights;
        storedHeights = heights;
        heights = tmp2;

        root = storedRoot;
    }

    /**
     * accept the stored state
     */
    @Override
    protected void acceptState() {
    } // nothing to do

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return the ith taxon in the list.
     */
    @Override
    public Taxon getTaxon(int taxonIndex) {
        return ((Node) getExternalNode(taxonIndex)).taxon;
    }

    // **************************************************************
    // MutableTaxonList IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    // ***********************************************************************
    // Private methods
    // ***********************************************************************

    // **************************************************************
    // Private inner classes
    // **************************************************************

    private class Node implements NodeRef {
        private Node(int number) {
            this(number, null);
        }

        private Node(int number, Taxon taxon) {
            this.number = number;
            this.taxon = taxon;
        }

        @Override
        public int getNumber() {
            return number;
        }

        @Override
        public void setNumber(int n) {
            throw new UnsupportedOperationException("Node is immutable");
        }

        public Taxon getTaxon() {
            return taxon;
        }

        final Taxon taxon;
        final int number;
    }
    // ***********************************************************************
    // Private members
    // ***********************************************************************


    /**
     * root node
     */
    private int root = -1;
    private int storedRoot;

    /**
     * list of internal nodes (including root)
     */
    private int[] edges = null;
    private int[] storedEdges = null;

    private double[] heights = null;
    private double[] storedHeights = null;

    private final NodeRef[] nodes;

    /**
     * number of nodes (including root and tips)
     */
    private final int nodeCount;

    /**
     * number of external nodes
     */
    private final int externalNodeCount;

    /**
     * number of internal nodes (including root)
     */
    private final int internalNodeCount;

    private boolean hasRates = false;
    private boolean hasTraits = false;

    @Override
    public boolean hasNodeTraits() {
        // TODO Auto-generated method stub
        return hasTraits;
    }

    @Override
    public void swapAllTraits(NodeRef node1, NodeRef node2) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("BFT model can not have traits! No trait to swap.");
    }

    @Override
    public boolean hasRates() {
        // TODO Auto-generated method stub
       return hasRates;
    }



}
