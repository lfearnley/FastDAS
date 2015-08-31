package gui;

import javax.swing.*;
import java.awt.event.*;
import java.io.File;

public class AutocuratorDialog extends JDialog {
	private JPanel contentPane;
	private JButton buttonOK;
	private JButton buttonCancel;
	private JTextField bioPaxFileField;
	private JButton selectBioPAXLevel3Button;
	private JTextField pathwayFileField;
	private JTextField outputDirectoryField;
	private JButton selectPathwayExclusionFileButton;
	private JButton selectOutputDirectoryButton;
	private JTextField speciesFileField;
	private JButton selectSpeciesExclusionFileButton;
	private JFileChooser fc = new JFileChooser();
	private File biopaxFile = null;
	private File pathwayFile = null;
	private File speciesFile = null;
	private File outputDirectory = null;

	public AutocuratorDialog() {
		setContentPane(contentPane);
		setModal(true);
		getRootPane().setDefaultButton(buttonOK);

		buttonOK.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onOK();
			}
		});

		buttonCancel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onCancel();
			}
		});

// call onCancel() when cross is clicked
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				onCancel();
			}
		});

// call onCancel() on ESCAPE
		contentPane.registerKeyboardAction(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onCancel();
			}
		}, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		selectBioPAXLevel3Button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				int selection = fc.showOpenDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					biopaxFile = fc.getSelectedFile();
					bioPaxFileField.setText(biopaxFile.getAbsolutePath());
				}
			}
		});
		selectPathwayExclusionFileButton.addActionListener(new ActionListener() {
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
		selectSpeciesExclusionFileButton.addActionListener(new ActionListener() {
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
		selectOutputDirectoryButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int selection = fc.showSaveDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					outputDirectory = fc.getSelectedFile();
					outputDirectoryField.setText(outputDirectory.getAbsolutePath());
				}
			}
		});
	}

	private void onOK() {
		if (this.biopaxFile == null || this.outputDirectory == null || this.speciesFile == null || this.pathwayFile ==
				null) {
			JOptionPane.showMessageDialog(null, "You must specify model, pathway and species exclusion files and an " +
					"output directory.");
		} else {
			FastDAS.MainClass.autoCuratorStart(biopaxFile, pathwayFile, speciesFile, outputDirectory);
			dispose();
		}
	}

	private void onCancel() {
		setVisible(false);
		dispose();
	}

	public static void main(String[] args) {
		AutocuratorDialog dialog = new AutocuratorDialog();
		dialog.pack();
		dialog.setVisible(true);
	}
}
