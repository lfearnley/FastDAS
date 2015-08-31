package gui;

import IO.BioPAXInputHandler;
import IO.TextOutputHandler;
import org.biopax.paxtools.model.Model;

import javax.swing.*;
import java.awt.event.*;
import java.io.File;

public class SpeciesGeneratorDialog extends JDialog {
	private JPanel contentPane;
	private JButton buttonOK;
	private JButton buttonCancel;
	private JTextField bioPaxFileTextField;
	private JTextField pathwayFileTextField;
	private JTextField outputDirectoryTextField;
	private JButton selectBioPAXLevel3Button;
	private JButton selectPathwayExclusionFileButton;
	private JButton selectOutputDirectoryButton;
	private JFileChooser fc = new JFileChooser();
	private File biopaxFile = null;
	private File pathwayFile = null;
	private File outputDirectory = null;

	public SpeciesGeneratorDialog() {
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
					bioPaxFileTextField.setText(biopaxFile.getAbsolutePath());
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
					pathwayFileTextField.setText(pathwayFile.getAbsolutePath());
				}
			}
		});
		selectOutputDirectoryButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int selection = fc.showOpenDialog(null);
				if (selection == JFileChooser.APPROVE_OPTION) {
					outputDirectory = fc.getSelectedFile();
					outputDirectoryTextField.setText(outputDirectory.getAbsolutePath());
				}
			}
		});
	}

	private void onOK() {
// add your code here
		if (outputDirectory == null || biopaxFile == null || pathwayFile == null) {
			JOptionPane.showMessageDialog(null, "You must specify model and pathway exclusion files and an output " +
					"directory.");
		} else {
			Model model = BioPAXInputHandler.readModelFromFile(biopaxFile);
			BioPAXInputHandler.removePathways(pathwayFile, model);
			TextOutputHandler.writeSaveFileForEditing(model, outputDirectory.getAbsolutePath() + System.getProperty
					("file.separator") + "species.txt");
			setVisible(false);
			dispose();
		}
	}

	private void onCancel() {
// add your code here if necessary
		setVisible(false);
		dispose();
	}

	public static void main(String[] args) {
		SpeciesGeneratorDialog dialog = new SpeciesGeneratorDialog();
		dialog.pack();
		dialog.setVisible(true);
	}
}
