package gui;

import FastDAS.KeyParam;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Created with IntelliJ IDEA.
 * User: uqlfearn
 * Date: 23/08/13
 * Time: 3:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class MainFrame {

	public static void main(String[] args) {
		//Schedule a job for the event-dispatching thread:
		//creating and showing this application's GUI.
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				createAndShowGUI();
			}
		});
	}

	private static void createAndShowGUI() {
		JFrame frame = new JFrame("LogicSim");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setJMenuBar(MainFrame.getMenuBar());
		frame.setContentPane(MainFrame.getContents());
		//Display the window.
		frame.setSize(450, 260);
		frame.setVisible(true);
	}

	private static Container getContents() {
		JPanel contentPane = new JPanel(new BorderLayout());
		contentPane.setOpaque(true);
		return contentPane;
	}

	private static JMenuBar getMenuBar() {
		JMenuBar menuBar = new JMenuBar();
		//Create the File menu and its menu items:
		JMenu fileMenu = new JMenu("File");
		menuBar.add(fileMenu);
		JMenuItem analyseMenuItem = new JMenuItem("Run Qualitative Logic System");
		fileMenu.add(analyseMenuItem);
		analyseMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				maininterface.main(new String[0]);
			}
		});
		JMenuItem autocuratorMenuItem = new JMenuItem("AutoCurator");
		fileMenu.add(autocuratorMenuItem);
		autocuratorMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				AutocuratorDialog.main(new String[0]);
			}
		});
		JMenuItem exitMenuItem = new JMenuItem("Exit");
		fileMenu.add(exitMenuItem);
		exitMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				System.exit(0);
			}
		});
		JMenu generateMenu = new JMenu("Generate");
		menuBar.add(generateMenu);
		JMenuItem pathwayListMenuItem = new JMenuItem("Generate Pathway List");
		generateMenu.add(pathwayListMenuItem);
		pathwayListMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				PathwayGeneratorDialog.main(new String[0]);
			}
		});
		JMenuItem speciesListMenuItem = new JMenuItem("Generate Species List");
		generateMenu.add(speciesListMenuItem);
		speciesListMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SpeciesGeneratorDialog.main(new String[0]);
			}
		});
		//Add the Menu Items to the menu:
		//Create the Help menu and its menu items:
		JMenu helpMenu = new JMenu("Help");
		menuBar.add(helpMenu);
		JMenuItem aboutMenuItem = new JMenuItem("About LogicSim");
		helpMenu.add(aboutMenuItem);
		helpMenu.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JOptionPane.showMessageDialog(new JFrame(), "LogicSim" + KeyParam.NEWLINE + "Version : " + KeyParam
						.VERSION_NUMBER + KeyParam.NEWLINE + "By Liam Fearnley");
			}
		});
		return menuBar;
	}

}
