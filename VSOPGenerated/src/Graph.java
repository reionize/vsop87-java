import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;
import java.awt.*;
import javax.swing.*;

public class Graph extends JFrame {
	public Graph() {
    JSplitPane sp = new JSplitPane();
    // not middle, just hasty setting
    sp.setDividerLocation(1228);
    add(sp);
    sp.setLeftComponent(new CirclePanel());
    sp.setRightComponent(new JPanel());
    setSize(1228 + 200,600);
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    setTitle("DEMO");
    setVisible(true);
    }
	
	public static void main(String[] args)  {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new Graph();
      }});
  }
	// inner class for demo, best not to use inner
	// for your purpose
	private class CirclePanel extends JPanel  {
		public CirclePanel() {
			setBackground(Color.WHITE);
		}
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
      // smoother look
			Graphics2D g2 = (Graphics2D) g;
			g2.setColor(Color.BLUE);
			double a = 2.6;
			double x;
			String[] l;
			try {
				Scanner s = new Scanner(new File("diff.txt"));
				while (s.hasNextLine()) {
					l = s.nextLine().split("    ");
					x = Double.parseDouble(l[0]);
					a = Double.parseDouble(l[1]);
					g2.fillOval((int)(2*(x - 1583)),200 - (int) ((a-2)*200),5,5);
					if (x % 2 == 0) {
						double y = 0.242/(1 + Math.pow(Math.E, -(x - 1717.689)/36.05)) + 2.450;
//						double y = 2.345322126959093E-4 * x + 2.1418026133261123;
//						double y = 0.5 * Math.pow((x - 1470.001)/2212.97, 0.45) + 2.4;
						g2.setColor(Color.RED);
						g2.fillOval((int)(2*(x - 1583)),200 - (int) ((y-2)*200),5,5);
						g2.setColor(Color.BLUE);
					}
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
//			g2.fillOval(500,600 - (int) (a*200),5,5);
		}
	}
}
