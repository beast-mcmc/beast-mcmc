package dr.evomodel.bigfasttree.thorney;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.xml.*;


// see NodeHeightProxyParameter.java 

public class RootHeightProxyParameter extends Parameter.Proxy implements Bounds<Double> {
    private TreeModel tree;
    public RootHeightProxyParameter(String name,
                                    TreeModel tree) {
        super(name, 1);
        this.tree = tree;
    }

    public TreeModel getTree() {
        return tree;
    }



    @Override
    public double getParameterValue(int dim) {
        return tree.getNodeHeight(tree.getRoot());
    }

    @Override
    public void setParameterValue(int dim, double value) {
        tree.setNodeHeight(tree.getRoot(), value);
        tree.pushTreeChangedEvent(tree.getRoot());
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        tree.setNodeHeightQuietly(tree.getRoot(), value);
    }
    @Override
    public Bounds<Double> getBounds() {
        return this;
    }
    //bounds implementation
    @Override
    public Double getUpperLimit(int dimension) {
            return Double.POSITIVE_INFINITY;
    }

    @Override
    public Double getLowerLimit(int dimension) {
        // TODO Auto-generated method stub
        NodeRef node = tree.getRoot();
            double max = 0.0;
            for(int i = 0; i < tree.getChildCount(node); i++) {
                max = Math.max(max, tree.getNodeHeight(tree.getChild(node,i)));
            }
            return max;
    }

    @Override
    public int getBoundsDimension() {
                return this.getDimension();       
    }
    public String toString() {
        StringBuilder buffer = new StringBuilder(String.valueOf(getParameterValue(0)));
        Bounds bounds = null;

        for (int i = 1; i < getDimension(); i++) {
            buffer.append("\t").append(String.valueOf(getParameterValue(i)));
        }
        return buffer.toString();
    }

    @Override
    public void fireParameterChangedEvent() {
        tree.pushTreeChangedEvent(TreeChangedEvent.create(true, true));
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    private static final String ROOT_HEIGHT_PARAMETER = "rootHeightProxyParameter";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            return new RootHeightProxyParameter(ROOT_HEIGHT_PARAMETER, tree);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return RootHeightProxyParameter.class;
        }

        @Override
        public String getParserName() {
            return ROOT_HEIGHT_PARAMETER;
        }
    };
}


