/*
 * RateCovarianceStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Statistic;
import dr.stats.DiscreteStatistics;

/**
 * A statistic that tracks the covariance of rates on branches
 *
 * @author Alexei Drummond
 */
public class RateCovarianceStatistic extends TreeStatistic {

    public RateCovarianceStatistic(String name, Tree tree, BranchRateModel branchRateModel) {
        super(name);
        this.tree = tree;
        this.branchRateModel = branchRateModel;

        int n = tree.getExternalNodeCount();
        childRate = new double[2 * n - 4];
        parentRate = new double[childRate.length];
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {

        int n = tree.getNodeCount();
        int index = 0;
        for (int i = 0; i < n; i++) {
            NodeRef child = tree.getNode(i);
            NodeRef parent = tree.getParent(child);
            if (parent != null && !tree.isRoot(parent)) {
                childRate[index] = branchRateModel.getBranchRate(tree, child);
                parentRate[index] = branchRateModel.getBranchRate(tree, parent);
                index += 1;
            }
        }
        return DiscreteStatistics.covariance(childRate, parentRate);
    }

    private Tree tree = null;
    private BranchRateModel branchRateModel = null;
    private double[] childRate = null;
    private double[] parentRate = null;
}
