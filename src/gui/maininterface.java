package gui;

import FastDAS.KeyParam;
import FastDAS.MainClass;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 15/08/13
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class maininterface {

	private static JFrame frame = new JFrame("FastDAS");

	public static void main(String[] args) {
		frame.setContentPane(new maininterface().mainform);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private JTextField bpFileField;
	private JButton bpL3BrowseButton;
	private JTextField pathwayFileField;
	private JButton pathwayBrowseButton;
	private JTextField speciesFileField;
	private JButton speciesBrowseButton;
	private JComboBox objectiveBox;
	private JComboBox optModeComboBox;
	private JComboBox outputTypeComboBox;
	private JButton RUNButton;
	private JButton cancelButton;
	private JTextField outputDirectoryField;
	private JButton outputDirectoryBrowseButton;
	private JPanel mainform;
	private JCheckBox saveMidConfidencePredictionsCheckBox;
	private JCheckBox saveLowConfidencePredictionsCheckBox;
	private JTextField textField1;
	private JTextField experimentDescriptorField;

	private File biopaxFile = null;
	private File speciesFile = null;
	private File pathwayFile = null;
	private File outputFile = null;
	private JFileChooser fc = new JFileChooser();

	public maininterface() {
		pathwayFileField.setEditable(false);
		bpFileField.setEditable(false);
		speciesFileField.setEditable(false);
		outputDirectoryField.setEditable(false);
		//Action listeners for the browse buttons in the MainClass.MainClass Interface Form
		bpL3BrowseButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				int selection = fc.showOpenDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					biopaxFile = fc.getSelectedFile();
					bpFileField.setText(biopaxFile.getAbsolutePath());
				}
			}
		});
		pathwayBrowseButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				int selection = fc.showOpenDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					pathwayFile = fc.getSelectedFile();
					pathwayFileField.setText(pathwayFile.getAbsolutePath());
				}
			}
		});
		speciesBrowseButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				int selection = fc.showOpenDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					speciesFile = fc.getSelectedFile();
					speciesFileField.setText(speciesFile.getAbsolutePath());
				}
			}
		});
		outputDirectoryBrowseButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int selection = fc.showSaveDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					outputFile = fc.getSelectedFile();
					if (outputFile.isDirectory()) {
						outputDirectoryField.setText(outputFile.getAbsolutePath());
					}
				}
			}
		});
		//Listener for cancel button - call System.exit for now, update with more graceful exit later.
		cancelButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				frame.setVisible(false); //Invis the JFrame
				frame.dispose(); //Dispose of the JFrame
			}
		});
		//Listener for OK button
		RUNButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				//To change body of implemented methods use File | Settings | File Templates.
				MainClass.dispatch(biopaxFile, pathwayFile, speciesFile, objectiveBox.getSelectedIndex(),
						optModeComboBox.getSelectedIndex(), outputTypeComboBox.getSelectedIndex(), outputFile,
						saveLowConfidencePredictionsCheckBox.isSelected(), saveMidConfidencePredictionsCheckBox
								.isSelected(), textField1.getText());
			}
		});
	}

	private void createUIComponents() {
		// TODO: place custom component creation code here
		objectiveBox = new JComboBox(KeyParam.objectiveOptions);
		optModeComboBox = new JComboBox(KeyParam.modeOptions);
		outputTypeComboBox = new JComboBox(KeyParam.saveOptions);
	}
}
