/*
 * ConstrainedTreeModel.java
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

package dr.evomodel.bigfasttree.thorney;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * A tree model that only allows node height operations an/or operations that maintain constrained clades
 * The tree structure is borrowed (copied) from a BigFastTree model.
 * TODO Could be refactored into an abstract class.
 */
public class ConstrainedTreeModel extends TreeModel {
    //
    // Public stuff
    //

    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";

    public ConstrainedTreeModel(Tree tree, Tree constraintsTree) {
        this(CONSTRAINED_TREE_MODEL, tree, false, false, constraintsTree);
    }

    public ConstrainedTreeModel(ConstrainedTreeModel baseTree){
        super(CONSTRAINED_TREE_MODEL,baseTree.isVariable());
        internalNodeCount = baseTree.getInternalNodeCount();
        externalNodeCount = baseTree.getExternalNodeCount();

        nodeCount = internalNodeCount + externalNodeCount;

        nodes = new NodeRef[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            NodeRef node = baseTree.getNode(i);
            if (baseTree.isExternal(node)) {
                nodes[i] = new Node(i, baseTree.getNodeTaxon(node));
            }else{
                nodes[i]= new Node(i);
            }
        }

        root = baseTree.getRoot().getNumber();

        edges = new int[nodeCount * 3];
        System.arraycopy(baseTree.edges,0,edges,0,nodeCount*3);

        storedEdges = new int[nodeCount * 3];

        heights = new double[nodeCount];
        System.arraycopy(baseTree.heights,0,heights,0,nodeCount);

        storedHeights = new double[nodeCount];

        //Now copy subtrees
        subtrees = new ArrayList<>();
        for (int i = 0; i < baseTree.getSubtreeCount(); i++) {
            WrappingSubtree baseSubtree = (WrappingSubtree) baseTree.getSubtree(i);
            HashSet<NodeRef> tips = new HashSet<>();
            for(int j=0; j<baseSubtree.getExternalNodeCount();j++){
                WrappingSubtree.Node baseSubtreeNode = (WrappingSubtree.Node) baseSubtree.getExternalNode(j);
                NodeRef myNode = nodes[baseSubtreeNode.getBaseNodeNumber()];
                tips.add(myNode);
            }
            WrappingSubtree subtree = new WrappingSubtree(this, tips, i);
            subtrees.add(subtree);
        }
    }

    public ConstrainedTreeModel(String name, Tree tree, Tree constraintsTree) {
        this(name, tree, false, false, constraintsTree);
    }


    public ConstrainedTreeModel(TreeModel tree, Tree constraintsTree) {
        this(CONSTRAINED_TREE_MODEL, tree, constraintsTree);
    }

    public ConstrainedTreeModel(String name, Tree tree, boolean fixHeights, boolean fixTree, Tree constraintsTree) {
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


        // need to make subtrees to nodes in the dataTree Subtrees
        Map<BitSet, NodeRef> constraintsTreeMap = getBitSetNodeMap(constraintsTree, constraintsTree);
        Map<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(constraintsTree, this);

        Map<NodeRef, NodeRef> constraintsNodeToTreeNode = new HashMap<>();
        for (BitSet clade :
                constraintsTreeMap.keySet()) {
            NodeRef constraintsTreeNode = constraintsTreeMap.get(clade);
            NodeRef treeModelNode = treeModelMap.get(clade);
            if(treeModelNode==null){
                throw new RuntimeException("All clades in the constraints tree must be present in the starting tree");
            }
            constraintsNodeToTreeNode.put(constraintsTreeNode, treeModelNode);
        }


        this.subtrees = new ArrayList<WrappingSubtree>();
        for (int i = 0; i < constraintsTree.getInternalNodeCount(); i++) {
            NodeRef constraintsTreeInternalNode = constraintsTree.getInternalNode(i);
            Set<NodeRef> tips = new HashSet<>();
            for (int j = 0; j < constraintsTree.getChildCount(constraintsTreeInternalNode); j++) {
                NodeRef child = constraintsTree.getChild(constraintsTreeInternalNode, j);
                tips.add(constraintsNodeToTreeNode.get(child));
            }
            WrappingSubtree subtree = new WrappingSubtree(this, tips, subtrees.size());
            subtrees.add(subtree);
        }

    }

    /**
     * Gets a HashMap of clade bitsets to nodes in tree. This is useful for comparing the topology of trees
     *
     * @param referenceTree the tree that will be used to define taxa and tip numbers
     * @param tree          the tree for which clades are being defined
     * @return A HashMap with a BitSet of descendent taxa as the key and a node as value
     */
    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree, Tree tree) {
        HashMap<BitSet, NodeRef> map = new HashMap<>();
        addBits(referenceTree, tree, tree.getRoot(), map);
        return map;
    }

    /**
     * A private recursive function used by getBitSetNodeMap
     * This is modeled after the addClades in CladeSet and getClades in compatibility statistic
     *
     * @param referenceTree the tree that will be used to define taxa and tip numbers
     * @param tree          the tree for which clades are being defined
     * @param node          current node
     * @param map           map that is being appended to
     */
    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap<BitSet, NodeRef> map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));

        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(addBits(referenceTree, tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }

    public int getSubtreeCount(){
        return subtrees.size();
    }

    public TreeModel getSubtree(NodeRef node,SubtreeContext context) {
        return subtrees.get(((Node) node).getSubtreeNumber(context));
    }


    public void copyEdgesAndHeights(ConstrainedTreeModel baseTree){
        // this assumes the two trees have the same subtrees
        System.arraycopy(baseTree.heights,0,heights,0,nodeCount);
        System.arraycopy(baseTree.edges,0,edges,0,nodeCount*3);


    }
    /**
     * This is a short cut function to return the subtree of node.
     * It defaults to the root context which returns the subtree that contains the node as a root if it is both
     * the root of one subtree and an external node in another. This context is almost exclusively used by other
     * classes.
     * @param node
     * @return
     */
    public TreeModel getSubtree(NodeRef node){
        return subtrees.get(((Node) node).getSubtreeNumber(SubtreeContext.IncludeRoot));
    }

    public TreeModel getSubtree(int i) {
        return subtrees.get(i);
    }

    private int getSubtreeIndex(NodeRef node,SubtreeContext context) {
        return ((Node) node).getSubtreeNumber(context);
    }

    public int getSubtreeIndex(NodeRef node) {
        return ((Node) node).getSubtreeNumber(SubtreeContext.IncludeRoot);
    }

    private NodeRef getNodeInSubtree(Tree tree, NodeRef nodeRef,SubtreeContext context) {
       return ((WrappingSubtree)tree).getWrappingNode(nodeRef,context);
    }
    public NodeRef getNodeInSubtree(Tree tree, NodeRef nodeRef) {
        return getNodeInSubtree( tree, nodeRef,SubtreeContext.IncludeRoot) ;
    }
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
     * tree.
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
        return 1.0;
        //mimics default tree model without node rates;
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
        if (isExternal(node)) {
            return 0;
        }
        int kids = 0;
        for (int i = 0; i < 2; i++) {
            if (edges[(node.getNumber() * 3) + i + 1] > -1) {
                kids += 1;
            }
            ;
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

    public void endTreeEditQuietly(){
        inEdit = false;
    }
    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    /**
     * Set a new node as root node.
     */
    @Override
    public void setRoot(NodeRef newRoot) {
        throw new RuntimeException("Cannot directly change the topology of a constrained tree! Please use a compatible operator");
    }
    public void setRootByForce(NodeRef node) {
        this.root = node.getNumber();
    }


    @Override
    public void addChild(NodeRef p, NodeRef c) {
        throw new RuntimeException("Cannot directly change the toplogy of a constrained tree! Please use a compatible operator");
    }

    private void addChildByForce(NodeRef p, NodeRef c) {
        addChildQuietlyByForce(p,c);
        pushTreeChangedEvent(TreeChangedEvent.create(p, false));
        pushTreeChangedEvent(TreeChangedEvent.create(c, false));
    }

    private void addChildQuietlyByForce(NodeRef p, NodeRef c){
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
    }


    @Override
    public void removeChild(NodeRef p, NodeRef c) {
        throw new RuntimeException("Cannot directly change the toplogy of a constrained tree! Please use a compatible operator");
    }

    private void removeChildQuietlyByForce(NodeRef p, NodeRef c){
        if (!inEdit) throw new RuntimeException("Must be in edit transaction to call this method!");

        int parent = p.getNumber();
        int child = c.getNumber();

        if (getChild(parent, 0) == child) {
            setChild(parent, 0, -1);
            //move other child up
            if (getChildCount(p) == 1) {
                setChild(parent, 0, getChild(parent, 1));
                setChild(parent, 1, -1);
            }

        } else if (getChild(parent, 1) == child) {
            setChild(parent, 1, -1);
        } else {
            throw new IllegalArgumentException("Child not in node");
        }
        setParent(child, -1);
    }
    private void removeChildByForce(NodeRef p, NodeRef c) {
        removeChildQuietlyByForce(p, c);
        pushTreeChangedEvent(TreeChangedEvent.create(p, false));
        pushTreeChangedEvent(TreeChangedEvent.create(c, false));
    }

    protected void clearTopology(){
        Arrays.fill(edges, -1);
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
     *
     * @param edges       Edges are provided as index: child number; parent: array entry
     * @param nodeHeights Also sets the node heights to the provided values
     * @param childOrder  Array that contains whether a child node is left or right child
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
                this.removeChildByForce(nodes[i], this.getChild(nodes[i], 0)); // nodes move into the first spot when it's filled
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
                    this.addChildByForce(this.getNode(edges[i]), this.getExternalNode(nodeMap[i]));
                    System.out.println("external: " + edges[i] + " > " + nodeMap[i]);
                } else {
                    //internal node
                    this.addChildByForce(this.getNode(edges[i]), this.getNode(i));
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
                if (i < this.externalNodeCount) {
                    if (childOrder[i] == 0 && getChild(nodes[edges[i]], 0) != nodes[nodeMap[i]]) {
                        //swap child nodes
                        NodeRef node = getNode(edges[i]);

                        NodeRef childOne = getChild(node, 0);
                        NodeRef childTwo = getChild(node, 1);

                        removeChildByForce(node, childOne);
                        removeChildByForce(node, childTwo);

                        addChildByForce(node, childTwo);
                        addChildByForce(node, childOne);
                    }
                } else {
                    if (childOrder[i] == 0 && getChild(nodes[edges[i]], 0) != nodes[i]) {
                        //swap child nodes
                        NodeRef node = getNode(edges[i]);

                        NodeRef childOne = getChild(node, 0);
                        NodeRef childTwo = getChild(node, 1);

                        removeChildByForce(node, childOne);
                        removeChildByForce(node, childTwo);

                        addChildByForce(node, childTwo);
                        addChildByForce(node, childOne);
                    }
                }
            }

        }

        this.setRootByForce(nodes[newRootIndex]);
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


    // ***********************************************************************
    // Private members
    // ***********************************************************************

    private final List<WrappingSubtree> subtrees;
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


    private static class Node implements NodeRef {
        private Node(int number) {
            this(number, null);
        }

        private Node(int number, Taxon taxon) {
            this.number = number;
            this.taxon = taxon;
            this.subtreeNumberIncludeRoots = -1;
            this.subtreeNumberIncludeTips = -1;
            this.subtreeNodeNumberIncludeRoots = -1;
            this.subtreeNodeNumberIncludeTips = -1;
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

        private int getSubtreeNumber(SubtreeContext context) {
            switch (context) {
                case IncludeRoot:
                    return subtreeNumberIncludeRoots;
                case IncludeTips:
                    return subtreeNumberIncludeTips;
                default:
                    throw new IllegalArgumentException("unknown subtree number context");
            }
        }

        private void setSubtreeNumber(int i, SubtreeContext context) {

            switch (context) {
                case IncludeRoot:
                    if (subtreeNumberIncludeRoots != -1) {
                        throw new UnsupportedOperationException("Node is immutable");
                    }
                    subtreeNumberIncludeRoots = i;
                    break;
                case IncludeTips:
                    if (subtreeNumberIncludeTips != -1) {
                        throw new UnsupportedOperationException("Node is immutable");
                    }
                    subtreeNumberIncludeTips = i;
                    break;
                default:
                    throw new IllegalArgumentException("unknown subtree number context");
            }
        }

        private int getSubtreeNodeNumber(SubtreeContext context) {
            switch (context) {
                case IncludeRoot:
                    return subtreeNodeNumberIncludeRoots;
                case IncludeTips:
                    return subtreeNodeNumberIncludeTips;
                default:
                    throw new IllegalArgumentException("unknown subtree number context");
            }
        }

        private void setSubtreeNodeNumber(int i, SubtreeContext context) {
            switch (context) {
                case IncludeRoot:
                    if (subtreeNodeNumberIncludeRoots != -1) {
                        throw new UnsupportedOperationException("Node is immutable");
                    }
                    subtreeNodeNumberIncludeRoots = i;
                    break;
                case IncludeTips:
                    if (subtreeNodeNumberIncludeTips != -1) {
                        throw new UnsupportedOperationException("Node is immutable");
                    }
                    subtreeNodeNumberIncludeTips = i;
                    break;
                default:
                    throw new IllegalArgumentException("unknown subtree number context");
            }
        }

        private final Taxon taxon;
        private final int number;
        private int subtreeNumberIncludeRoots;
        private int subtreeNumberIncludeTips;
        private int subtreeNodeNumberIncludeRoots;
        private int subtreeNodeNumberIncludeTips;

    }

    /**
     * A private enum for define in what context we are mapping a Constraints tree node to a subtree node.
     * The root of subtrees are also tips in their parent tree. The public api always used the Include root
     * context. In these cases a node in the base tree that is in two subtrees belongs to the tree where it is
     * the root.
     */
    protected enum SubtreeContext {
        IncludeRoot,
        IncludeTips
    }

    /**
     * A private tree model class that wraps the constrained tree. It does not have it's own child
     * parent mapping but instead relies on the constrained tree. The node mapping between node in one of these trees and
     * the constrained tree is constant and final. The topology of these subtrees can be operated on which changes
     * the topology in the constrained tree, without breaking any constraints.
     */
    protected class WrappingSubtree extends TreeModel {

        public static final String WRAPPED_TREE_MODEL = "wrappedTreeModel";

        public WrappingSubtree(String name, ConstrainedTreeModel tree, Set<NodeRef> externalNodes, int number) {
            super(name, tree.isVariable());
            setId(name);

            this.wrappedTree = tree;
            this.number = number;
            externalNodeCount = externalNodes.size();
            nodeCount = externalNodeCount * 2 - 1;
            internalNodeCount = nodeCount - externalNodeCount;

            nodes = new Node[nodeCount];
            int[] wrappedNodes = new int[externalNodeCount];

            int currentNodeIndex = 0;
            for (NodeRef wrappedNode : externalNodes) {
                nodes[currentNodeIndex] = new Node(currentNodeIndex, wrappedNode.getNumber());
                ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNumber(number,SubtreeContext.IncludeTips);
                ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNodeNumber(currentNodeIndex,SubtreeContext.IncludeTips);
                //if these are true tip then we will never see them again
                if(tree.isExternal(wrappedNode)){
                    ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNumber(number,SubtreeContext.IncludeRoot);
                    ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNodeNumber(currentNodeIndex,SubtreeContext.IncludeRoot);
                }
                wrappedNodes[currentNodeIndex] = wrappedNode.getNumber();
                currentNodeIndex++;
            }

            NodeRef wrappedNode = TreeUtils.getCommonAncestor(tree, wrappedNodes);

            Node rootNode = new Node(currentNodeIndex, wrappedNode.getNumber());
            nodes[currentNodeIndex]=rootNode;
            root=currentNodeIndex;
            ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNumber(number,SubtreeContext.IncludeRoot);
            ((ConstrainedTreeModel.Node) wrappedNode).setSubtreeNodeNumber(currentNodeIndex,SubtreeContext.IncludeRoot);
            currentNodeIndex++;
            for (int j = 0; j < tree.getChildCount(wrappedNode); j++) {
                NodeRef childToWrap = tree.getChild(wrappedNode, j);
                currentNodeIndex= traverseAndSetup(tree, childToWrap, externalNodes, currentNodeIndex);
            }


        }

        public WrappingSubtree(ConstrainedTreeModel tree, Set<NodeRef> externalNodes, int number) {
            this(WRAPPED_TREE_MODEL, tree, externalNodes,number);
        }


        private int traverseAndSetup(Tree tree, NodeRef node, Set<NodeRef> externalNodes, int currentIndex) {
            if (tree.isExternal(node)) {
                if (!externalNodes.contains(node)) {
                    throw new IllegalArgumentException("Wrapped paraphyletic wrapped trees are not supported");
                }
            }
            if (!externalNodes.contains(node)) {

                nodes[currentIndex] = new Node(currentIndex, node.getNumber());
                ((ConstrainedTreeModel.Node) node).setSubtreeNumber(number,SubtreeContext.IncludeTips);
                ((ConstrainedTreeModel.Node) node).setSubtreeNodeNumber(currentIndex,SubtreeContext.IncludeTips);
                ((ConstrainedTreeModel.Node) node).setSubtreeNumber(number,SubtreeContext.IncludeRoot);
                ((ConstrainedTreeModel.Node) node).setSubtreeNodeNumber(currentIndex,SubtreeContext.IncludeRoot);
                currentIndex++;
                //new counter to count decendents;
                for (int i = 0; i < tree.getChildCount(node); i++) {
                    currentIndex = traverseAndSetup(tree, tree.getChild(node, i), externalNodes, currentIndex);
                }
            }
            return currentIndex;
        }

        public NodeRef getNodeInWrappedTree(NodeRef node) {
            return wrappedTree.getNode(((Node) node).baseNodeNumber);
        }

        protected NodeRef getWrappingNode(NodeRef node, SubtreeContext context) {
            int index = ((ConstrainedTreeModel.Node) node).getSubtreeNodeNumber(context);
            return nodes[index];
        }

        /**
         * Add child to the children of parent.
         *
         * @param parent
         * @param child
         * @throws IllegalArgumentException If child is already a child of parent
         */
        @Override
        public void addChild(NodeRef parent, NodeRef child) {
            NodeRef wrappedParent = getNodeInWrappedTree(parent);
            NodeRef wrappedChild = getNodeInWrappedTree(child);

            if (isRoot(child)) {
                NodeRef rootParent = wrappedTree.getParent(wrappedChild);
                wrappedTree.removeChildByForce(rootParent, wrappedChild);
            }
            wrappedTree.addChildByForce(wrappedParent, wrappedChild);
        }

        public void addChildQuietly(NodeRef parent, NodeRef child) {
            NodeRef wrappedParent = getNodeInWrappedTree(parent);
            NodeRef wrappedChild = getNodeInWrappedTree(child);

            if (isRoot(child)) {
                NodeRef rootParent = wrappedTree.getParent(wrappedChild);
                wrappedTree.removeChildQuietlyByForce(rootParent, wrappedChild);
            }
            wrappedTree.addChildQuietlyByForce(wrappedParent, wrappedChild);
        }

        /**
         * Removes child from the children of parent.
         *
         * @param parent
         * @param child
         * @throws IllegalArgumentException If child is not a child of parent
         */
        @Override
        public void removeChild(NodeRef parent, NodeRef child) {
            NodeRef wrappedParent = getNodeInWrappedTree(parent);
            NodeRef wrappedChild = getNodeInWrappedTree(child);
            wrappedTree.removeChildByForce(wrappedParent, wrappedChild);
        }
        public void removeChildQuietly(NodeRef parent, NodeRef child) {
            NodeRef wrappedParent = getNodeInWrappedTree(parent);
            NodeRef wrappedChild = getNodeInWrappedTree(child);
            wrappedTree.removeChildQuietlyByForce(wrappedParent, wrappedChild);
        }
        /**
         * Replace child with another
         *
         * @param node
         * @param child    of node to replace
         * @param newChild replacment child
         */
        @Override
        public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
            NodeRef wrappedParent = getNodeInWrappedTree(node);
            NodeRef wrappedChild = getNodeInWrappedTree(child);
            NodeRef wrappedNewChild = getNodeInWrappedTree(newChild);
            wrappedTree.replaceChild(wrappedParent, wrappedChild, wrappedNewChild);
        }

        /**
         * set the height of the ith node in the tree (where the first n are internal).
         *
         * @param node
         * @param height
         */
        @Override
        public void setNodeHeight(NodeRef node, double height) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            wrappedTree.setNodeHeight(wrappedNode, height);
        }


        /**
         * set the rate of the ith node in the tree (where the first n are internal).
         *
         * @param node
         * @param rate
         */
        @Override
        public void setNodeRate(NodeRef node, double rate) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            wrappedTree.setNodeRate(wrappedNode, rate);
        }

        /**
         * Sets an named attribute for a given node.
         *
         * @param node  the node whose attribute is being set.
         * @param name  the name of the attribute.
         * @param value the new value of the attribute.
         */
        @Override
        public void setNodeAttribute(NodeRef node, String name, Object value) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            wrappedTree.setNodeAttribute(wrappedNode, name, value);
        }

        @Override
        public void setMultivariateTrait(NodeRef node, String name, double[] value) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            wrappedTree.setMultivariateTrait(wrappedNode, name, value);
        }

        /**
         * Modifies the current tree by adopting the provided collection of edges
         *
         * @param edges       Edges are provided as index: child number; parent: array entry
         * @param nodeHeights Also sets the node heights to the provided values
         * @param childOrder  Array that contains whether a child node is left or right child
         * @param taxaNames
         */
        @Override
        public void adoptTreeStructure(int[] edges, double[] nodeHeights, int[] childOrder, String[] taxaNames) {
//Nothing to do I think?
        }

        @Override
        public boolean isTreeValid() {
            return wrappedTree.isTreeValid();
        }

        public void setRoot(NodeRef node) {
            throw new UnsupportedOperationException("wrapped trees can not change roots, yet");
        }

        public NodeRef getRoot() {
            return getNode(root);
        }

        public boolean isRoot(NodeRef nodeRef) {
            return nodeRef.getNumber() == root;
        }


        @Override
        public void setNodeHeightQuietly(NodeRef node, double height) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            wrappedTree.setNodeHeightQuietly(wrappedNode, height);
        }


        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {

        }

        /**
         * This method is called whenever a parameter is changed.
         * <p/>
         * It is strongly recommended that the model component sets a "dirty" flag and does no
         * further calculations. Recalculation is typically done when the model component is asked for
         * some information that requires them. This mechanism is 'lazy' so that this method
         * can be safely called multiple times with minimal computational cost.
         *
         * @param variable
         * @param index
         * @param type
         */
        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        }

        /**
         * This call specifies that the current state is accept. Most models will not need to do anything.
         * Sub-models are handled automatically and do not need to be considered in this method.
         */
        @Override
        protected void acceptState() {

        }

        private final int root;


        /**
         * set the length of the branch above the ith node in the tree (where the first n are internal).
         *
         * @param node
         * @param length
         */
        @Override
        public void setBranchLength(NodeRef node, double length) {
            throw new UnsupportedOperationException("WrappedTreeModels cannot have branch lengths set");
        }


        @Override
        public double[] getMultivariateNodeTrait(NodeRef node, String name) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            return wrappedTree.getMultivariateNodeTrait(wrappedNode, name);
        }


        /**
         * @return a count of the number of nodes (internal + external) in this
         * tree, currently connected from the root node.
         */
        @Override
        public int getNodeCount() {
            return nodeCount;
        }

        /**
         * @param i node index, terminal nodes are first
         * @return the ith node.
         */
        @Override
        public NodeRef getNode(int i) {
            return nodes[i];
        }

        /**
         * @param i index of an internal node
         * @return the ith internal node.
         */
        @Override
        public NodeRef getInternalNode(int i) {
            return nodes[externalNodeCount + i];
        }

        /**
         * @param i the index of an external node
         * @return the ith external node.
         */
        @Override
        public NodeRef getExternalNode(int i) {
            return nodes[i];
        }

        /**
         * @return a count of the number of external nodes (tips) in this
         * tree, currently connected from the root node.
         */
        @Override
        public int getExternalNodeCount() {
            return externalNodeCount;
        }

        /**
         * @return a count of the number of internal nodes in this
         * tree, currently connected from the root node.
         */
        @Override
        public int getInternalNodeCount() {
            return internalNodeCount;
        }

        /**
         * @param node the node to retrieve the taxon of
         * @return the taxon of this node.
         */
        @Override
        public Taxon getNodeTaxon(NodeRef node) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            return wrappedTree.getNodeTaxon(wrappedNode);
        }

        /**
         * @param node the node to retrieve height of
         * @return the height of node in the tree.
         */
        @Override
        public double getNodeHeight(NodeRef node) {
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            return wrappedTree.getNodeHeight(wrappedNode);
        }

        /**
         * @param node the node to retrieve the rate of
         * @return the rate of node in the tree.
         */
        @Override
        public double getNodeRate(NodeRef node) {
            NodeRef wrappedNode = getNodeInWrappedTree((Node) node);
            return wrappedTree.getNodeRate(wrappedNode);
        }

        /**
         * @param node the node whose attribute is being fetched.
         * @param name the name of the attribute of interest.
         * @return an object representing the named attributed for the given node.
         */
        @Override
        public Object getNodeAttribute(NodeRef node, String name) {
            NodeRef wrappedNode = getNodeInWrappedTree((Node) node);
            return wrappedTree.getNodeAttribute(wrappedNode, name);
        }

        /**
         * @param node the node whose attribute is being fetched.
         * @return an interator of attribute names available for this node.
         */
        @Override
        public Iterator getNodeAttributeNames(NodeRef node) {
            NodeRef wrappedNode = getNodeInWrappedTree((Node) node);
            return wrappedTree.getNodeAttributeNames(wrappedNode);
        }


        /**
         * @param node the node to test if external
         * @return whether the node is external.
         */
        @Override
        public boolean isExternal(NodeRef node) {
            return node.getNumber() < externalNodeCount;
        }


        /**
         * @param node the node to get child count of
         * @return the number of children of node.
         */
        @Override
        public int getChildCount(NodeRef node) {
            if (isExternal(node)) {
                return 0;
            }
            return wrappedTree.getChildCount(getNodeInWrappedTree(node));
        }


        /**
         * @param node the node to get jth child of
         * @param j    the index of child to retrieve
         * @return the jth child of node
         */
        @Override
        public NodeRef getChild(NodeRef node, int j) {
            if (isExternal(node)) {
                throw new IllegalArgumentException("Attempted to access a child of an external node in a wrapped tree");
            }
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            NodeRef wrappedChild = wrappedTree.getChild(wrappedNode, j);

            return getWrappingNode(wrappedChild, SubtreeContext.IncludeTips);
        }

        @Override
        public NodeRef getParent(NodeRef node) {
            if (isRoot(node)) {
                return null;
            }
            NodeRef wrappedNode = getNodeInWrappedTree(node);
            NodeRef wrappedParent = wrappedTree.getParent(wrappedNode);
            return getWrappingNode(wrappedParent, SubtreeContext.IncludeRoot);
        }


        /**
         * @param taxonIndex
         * @return the ith taxon.
         */
        @Override
        public Taxon getTaxon(int taxonIndex) {
            int wrappedTaxonIndex = nodes[taxonIndex].baseNodeNumber;
            return wrappedTree.getTaxon(wrappedTaxonIndex);
        }

        @Override
        public NodeRef[] getNodes() {
            return nodes;
        }

        @Override
        public boolean beginTreeEdit() {
            wrappedTree.beginTreeEdit();
            return super.beginTreeEdit();
        }

        @Override
        public void endTreeEdit() {
            wrappedTree.endTreeEdit();
            super.endTreeEdit();
        }


        /**
         * Additional state information, outside of the sub-model is stored by this call.
         */
        @Override
        protected void storeState() {
        }


        /**
         * After this call the model is guaranteed to have returned its extra state information to
         * the values coinciding with the last storeState call.
         * Sub-models are handled automatically and do not need to be considered in this method.
         */
        @Override
        protected void restoreState() {
        }

        private final ConstrainedTreeModel wrappedTree;
        private final Node[] nodes;
        private final int nodeCount;
        private final int externalNodeCount;
        private final int internalNodeCount;
        private final int number;

        private class Node implements NodeRef {
            private Node(int number, int baseNodeNumber) {
                this.number = number;
                this.baseNodeNumber = baseNodeNumber;
            }

            @Override
            public int getNumber() {
                return number;
            }

            @Override
            public void setNumber(int n) {
                throw new UnsupportedOperationException("Node is immutable");
            }

            public int getBaseNodeNumber() {
                return baseNodeNumber;
            }

            public final int number;
            public final int baseNodeNumber;

        }

        @Override
        public boolean hasNodeTraits() {
            // TODO Auto-generated method stub
            throw new UnsupportedOperationException("Unimplemented method 'hasNodeTraits'");
        }

        @Override
        public void swapAllTraits(NodeRef node1, NodeRef node2) {
            // TODO Auto-generated method stub
            throw new UnsupportedOperationException("Unimplemented method 'swapAllTraits'");
        }

        @Override
        public boolean hasRates() {
            // TODO Auto-generated method stub
            throw new UnsupportedOperationException("Unimplemented method 'hasRates'");
        }
    }

    @Override
    public boolean hasNodeTraits() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'hasNodeTraits'");
    }

    @Override
    public void swapAllTraits(NodeRef node1, NodeRef node2) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'swapAllTraits'");
    }

    @Override
    public boolean hasRates() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'hasRates'");
    }

}
