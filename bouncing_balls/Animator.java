package bouncing_balls;

import java.awt.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.Timer;

import bouncing_balls.Model.Ball;


@SuppressWarnings("serial")
public final class Animator extends JPanel implements ActionListener {

	public Animator(int pixelWidth, int pixelHeight, int fps) {
		super(true);
		this.timer = new Timer(1000 / fps, this);
		this.deltaT = 1.0 / fps;
		this.model = new Model(pixelWidth / pixelsPerMeter, pixelHeight / pixelsPerMeter);
		this.setOpaque(false);
		this.setPreferredSize(new Dimension(pixelWidth, pixelHeight));
	}

	private static final double pixelsPerMeter = 200;
	private Model model;
	private Timer timer;
	private double deltaT;

	public void start() {
		timer.start();
	}

	public void stop() {
    	timer.stop();
    }

	@Override
	protected void paintComponent(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		// clear the canvas
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, this.getWidth(), this.getHeight());
		// draw balls
		g2.setColor(Color.RED);

		for (Ball b : model.balls) {
			double x = b.x - b.radius;
			double y = b.y + b.radius;
			// paint balls (y-coordinates are inverted)
			Ellipse2D.Double e = new Ellipse2D.Double(
					x * pixelsPerMeter,
					this.getHeight() - (y * pixelsPerMeter),
					b.radius * 2 * pixelsPerMeter,
					b.radius * 2 * pixelsPerMeter);
			g2.fill(e);
		}
		g2.setColor(Color.BLACK);
		//draw Rectangles
		for(Model.Rectangle r : model.rectangles){

			Path2D.Double p = new Path2D.Double();

			p.moveTo(r.x1 * pixelsPerMeter,this.getHeight() - (r.y1 * pixelsPerMeter));
			p.lineTo(r.x3 * pixelsPerMeter,this.getHeight() - (r.y3 * pixelsPerMeter));
			p.lineTo(r.x4 * pixelsPerMeter,this.getHeight() - (r.y4 * pixelsPerMeter));
			p.lineTo(r.x2 * pixelsPerMeter,this.getHeight() - (r.y2 * pixelsPerMeter));
			p.closePath();

			g2.setColor(Color.red);
			g2.fill(p);
		}
		Toolkit.getDefaultToolkit().sync();
	}

    @Override
    public void actionPerformed(ActionEvent e) {
    	model.step(deltaT);
    	this.repaint();
    }

	public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                Animator anim = new Animator(1920, 1080, 120);
                JFrame frame = new JFrame("Bouncing balls");
            	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            	frame.add(anim);
            	frame.pack();
            	frame.setLocationRelativeTo(null);
            	frame.setVisible(true);
            	anim.start();
            }
        });
    }
}
