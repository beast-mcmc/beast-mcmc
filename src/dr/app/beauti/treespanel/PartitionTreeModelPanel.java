/*
 * PartitionTreeModelPanel.java
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

package dr.app.beauti.treespanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.TreeHolder;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreeAsDataType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.BEAUTiImporter;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionTreeModelPanel extends OptionsPanel {

    private static final long serialVersionUID = 8096349200725353543L;

    private static final Insets FIELD_INSETS = new Insets(3,6,3,6);

    private final String NO_TREE = "no tree loaded";


    private final ButtonGroup startingTreeGroup = new ButtonGroup();
    private final JRadioButton randomTreeRadio = new JRadioButton("Random starting tree");
    private final JRadioButton upgmaTreeRadio = new JRadioButton("UPGMA starting tree");
    private final JRadioButton userTreeRadio = new JRadioButton("User-specified starting tree");
    private final JLabel userTreeLabel = new JLabel("User-specified tree:");
    private final JComboBox userTreeCombo = new JComboBox();
    private final JLabel userTreeInfo = new JLabel("<html>" +
            "Use a tree imported using the 'Import Data' menu option.<br>" +
            "Starting trees that are not rooted and strictly bifurcating (binary) will be randomly resolved.</html>");


    private JComboBox treeAsDataModelCombo = new JComboBox();
    private RealNumberField thorneyGenomeLength = new RealNumberField(1.0, Double.MAX_VALUE);

    private final OptionsPanel thorneyBEASTPanel;
    private final JLabel thorneyBEASTInfo = new JLabel("<html>" +
            "Use the tree and branch lengths in substitutions per site as data, integrating over unresolved<br> " +
            "polytomies and sampling branch lengths.<br><br></html>");

    private final JLabel thorneyBEASTCitation = new JLabel("<html>"+
    "Citation: du Plessis L, McCrone JT, <i>et al.</i> (2021) Establishment and lineage dynamics of the SARS-CoV-2 epidemic<br>" +
    "in the UK. <i>Science</i> <b>371</b>, 708-712. DOI:10.1126/science.abf2946 </html>");

    private final OptionsPanel empiricalTreePanel;
    private final JLabel empiricalTreeLabel = new JLabel("Trees filename:");
    private final JCheckBox empiricalExternalFileCheck = new JCheckBox("Read empirical trees from an external file:");
    private final JTextArea empiricalFilenameField = new JTextArea("empirical.trees");

    private final JLabel empiricalTreeInfo = new JLabel(
            "<html>" +
                    "Specify a file from which to load the empirical trees when BEAST is run. The trees should<br>" +
                    "have the same taxa as the loaded tree partition.<br>" +
                    "BEAST will look for the file in the same folder as the XML file or in the current working<br>" +
                    "directory. The file can also be specified with a relative or absolute path." +
                    "</html>");

    private final RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.POSITIVE_INFINITY, "Init root height");

    private BeautiOptions options = null;
    private final BeautiFrame parent;
    private boolean settingOptions = false;

    PartitionTreeModel partitionTreeModel;

    public PartitionTreeModelPanel(final BeautiFrame parent, PartitionTreeModel partitionTreeModel, final BeautiOptions options) {
        super(12, (OSType.isMac() ? 6 : 24));

        this.partitionTreeModel = partitionTreeModel;
        this.options = options;
        this.parent = parent;

        PanelUtils.setupComponent(initRootHeightField);
        initRootHeightField.setColumns(10);
        initRootHeightField.setEnabled(false);

        PanelUtils.setupComponent(randomTreeRadio);
        PanelUtils.setupComponent(upgmaTreeRadio);
        PanelUtils.setupComponent(userTreeRadio);

        startingTreeGroup.add(randomTreeRadio);
        startingTreeGroup.add(upgmaTreeRadio);
        startingTreeGroup.add(userTreeRadio);

        randomTreeRadio.setSelected(this.partitionTreeModel.getStartingTreeType() == StartingTreeType.RANDOM);
        upgmaTreeRadio.setSelected(this.partitionTreeModel.getStartingTreeType() == StartingTreeType.UPGMA);
        userTreeRadio.setSelected(this.partitionTreeModel.getStartingTreeType() == StartingTreeType.USER);
        userTreeRadio.setEnabled(!options.userTrees.isEmpty());

        boolean enabled = this.partitionTreeModel.getStartingTreeType() == StartingTreeType.USER;
        userTreeLabel.setEnabled(enabled);
        userTreeCombo.setEnabled(enabled);
        userTreeInfo.setEnabled(enabled);

       for (TreeAsDataType treeAsDataType : TreeAsDataType.values()) {
           treeAsDataModelCombo.addItem(treeAsDataType);
       }
        PanelUtils.setupComponent(treeAsDataModelCombo);
        treeAsDataModelCombo.addItemListener(ev -> {
            this.partitionTreeModel.setTreeAsDataType((TreeAsDataType)treeAsDataModelCombo.getSelectedItem());
            setupPanel();
        });

        ActionListener listener = actionEvent -> {
            if (randomTreeRadio.isSelected()) {
                this.partitionTreeModel.setStartingTreeType(StartingTreeType.RANDOM);
            } else if (upgmaTreeRadio.isSelected()) {
                this.partitionTreeModel.setStartingTreeType(StartingTreeType.UPGMA);
            } else if (userTreeRadio.isSelected()) {
                this.partitionTreeModel.setStartingTreeType(StartingTreeType.USER);
            }
            boolean enabled1 = this.partitionTreeModel.getStartingTreeType() == StartingTreeType.USER;
            userTreeLabel.setEnabled(enabled1);
            userTreeCombo.setEnabled(enabled1);
            userTreeInfo.setEnabled(enabled1);

            parent.setDirty();
        };
        randomTreeRadio.addActionListener(listener);
        upgmaTreeRadio.addActionListener(listener);
        userTreeRadio.addActionListener(listener);

        PanelUtils.setupComponent(userTreeCombo);

        userTreeCombo.addItemListener(ev -> setUserSpecifiedStartingTree());

        thorneyBEASTPanel = new OptionsPanel();
        thorneyBEASTPanel.setOpaque(false);

        PanelUtils.setupComponent(thorneyGenomeLength);

        thorneyGenomeLength.setToolTipText("<html>" +
                "This scales branches in the input tree <br>" +
                "to represent expected number of mutations in the genome.<br>" +
                "It should likely be the length of the alignment to used to build the tree.<html>");
        thorneyGenomeLength.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                PartitionTreeModelPanel.this.partitionTreeModel.setThorneyScaler(thorneyGenomeLength.getValue());
                parent.setDirty();
            }
        });
        thorneyGenomeLength.setColumns(10);


        PanelUtils.setupComponent(empiricalExternalFileCheck);
        empiricalExternalFileCheck.addActionListener(ev -> {
            empiricalTreeInfo.setEnabled(empiricalExternalFileCheck.isSelected());
            empiricalTreeLabel.setEnabled(empiricalExternalFileCheck.isSelected());
            empiricalFilenameField.setEnabled(empiricalExternalFileCheck.isSelected());
            this.partitionTreeModel.setUsingExternalEmpiricalTreeFile(empiricalExternalFileCheck.isSelected());
            parent.setDirty();
        });

        empiricalFilenameField.setBorder(new EmptyBorder(FIELD_INSETS));
        empiricalTreePanel = new OptionsPanel();
        empiricalTreePanel.setOpaque(false);
        empiricalFilenameField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                PartitionTreeModelPanel.this.partitionTreeModel.setEmpiricalTreesFilename(empiricalFilenameField.getText().trim());
                parent.setDirty();
            }
        });

        setupPanel();
        setOptions();
    }

    private void setUserSpecifiedStartingTree() {
        if (userTreeCombo.getSelectedItem() != null && (!userTreeCombo.getSelectedItem().toString().equalsIgnoreCase(NO_TREE))) {
            Tree selectedTree = getSelectedUserTree();
            if (selectedTree != null) {
                partitionTreeModel.setUserStartingTree(selectedTree);
            } else {
                JOptionPane.showMessageDialog(parent, "The selected user-specified starting tree " +
                                "is not fully bifurcating.\nBEAST requires rooted, bifurcating (binary) trees.",
                        "Illegal user-specified starting tree",
                        JOptionPane.ERROR_MESSAGE);

                userTreeCombo.setSelectedItem(NO_TREE);
                partitionTreeModel.setUserStartingTree(null);
            }
        }
    }

    public void setupPanel() {

        removeAll();

        if (partitionTreeModel.getDataType().getType() != DataType.TREE) {
            addSpanningComponent(randomTreeRadio);
            addSpanningComponent(upgmaTreeRadio);
            addSpanningComponent(userTreeRadio);

            addComponents(userTreeLabel, userTreeCombo);
            userTreeCombo.removeAllItems();
            if (options.userTrees.isEmpty()) {
                userTreeCombo.addItem(NO_TREE);
            } else {
                Object selectedItem = userTreeCombo.getSelectedItem();
                for (TreeHolder tree : options.userTrees.values()) {
                    userTreeCombo.addItem(tree);
                }
                if (selectedItem != null) {
                    userTreeCombo.setSelectedItem(selectedItem);
                } else {
                    userTreeCombo.setSelectedIndex(0);
                }
            }

            addComponent(userTreeInfo);

        } else {

            addComponentWithLabel("Tree as data model:", treeAsDataModelCombo);

            thorneyBEASTPanel.removeAll();
            thorneyBEASTPanel.addComponent(thorneyBEASTInfo);
            thorneyBEASTPanel.addComponentWithLabel("Scaler", thorneyGenomeLength);
            thorneyBEASTPanel.addComponent(thorneyBEASTCitation);

            empiricalTreePanel.removeAll();
            empiricalTreePanel.addComponent(new JLabel( "Use the loaded trees as an empirical set to sample over."));
            TreeHolder trees = partitionTreeModel.getTreePartitionData().getTrees();
            empiricalTreePanel.addComponent(new JLabel(""));
            empiricalTreePanel.addComponent(new JLabel("<html>" +
                    "This tree partition, " + partitionTreeModel.getName() + ", is a set of " + trees.getTreeCount() + " trees.<br>" +
                    "<br>" +
                    (trees.getTrees().size() == 1 ? "Only one tree loaded - specify file to read full tree set from:": "") +
                    "</html>"));

            empiricalTreePanel.addComponent(empiricalExternalFileCheck);
            empiricalExternalFileCheck.setSelected(trees.getTrees().size() == 1);

            empiricalTreeInfo.setEnabled(empiricalExternalFileCheck.isSelected());
            empiricalTreeLabel.setEnabled(empiricalExternalFileCheck.isSelected());
            empiricalFilenameField.setEnabled(empiricalExternalFileCheck.isSelected());

            empiricalTreePanel.addComponent(empiricalTreeInfo);
            empiricalTreePanel.addComponents(empiricalTreeLabel, empiricalFilenameField);
            empiricalFilenameField.setColumns(32);
            empiricalFilenameField.setEditable(true);

            switch ((TreeAsDataType)treeAsDataModelCombo.getSelectedItem()) {
                case EMPRICAL_TREES:
                    addSpanningComponent(empiricalTreePanel);
                    break;
                case THORNEY_BEAST:
                    addSpanningComponent(thorneyBEASTPanel);
                    break;
            }
        }

        validate();
        repaint();
    }

    public void setOptions() {

        if (partitionTreeModel == null) {
            return;
        }

        settingOptions = true;
        initRootHeightField.setValue(partitionTreeModel.getInitialRootHeight());

        userTreeRadio.setEnabled(!options.userTrees.isEmpty());

        settingOptions = false;

        treeAsDataModelCombo.setSelectedItem(partitionTreeModel.getTreeAsDataType());
        empiricalExternalFileCheck.setSelected(partitionTreeModel.isUsingExternalEmpiricalTreeFile());
        empiricalFilenameField.setText(partitionTreeModel.getEmpiricalTreesFilename());

        
        thorneyGenomeLength.setValue(partitionTreeModel.getThorneyScaler());
        

    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;
    }

    public boolean isBifurcatingTree(Tree tree, NodeRef node) {
        if (tree.getChildCount(node) > 2) return false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (!isBifurcatingTree(tree, tree.getChild(node, i))) return false;
        }
        return true;
    }

    private Tree getSelectedUserTree() {
        TreeHolder treeHolder = (TreeHolder) userTreeCombo.getSelectedItem();
        return treeHolder.getTrees().get(0);
    }

}