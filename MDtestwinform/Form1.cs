using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MDtestwinform
{
    public partial class Form1 : Form
    {
        Thread listen;
        delegate void Change(int progress, List<Node> nodes, VecR region, int step);
        public Form1()
        {
            InitializeComponent();
        }

        
        private void Update(int per, List<Node> nodes, VecR region,int step)
        {
            progressBar1.Value = per;
            this.Text = per.ToString();
            Bitmap bmp = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            Graphics g = Graphics.FromImage(bmp);
            int radius = 10;
            foreach(Node node in nodes)
            {
                float x = (float)((node.Getr.x + 0.5 * region.x) * pictureBox1.Width / region.x);
                float y = (float)((node.Getr.y + 0.5 * region.y) * pictureBox1.Height / region.y);
                g.DrawEllipse(new Pen(Color.Black), x, y, radius, radius);
          //      Console.WriteLine(x);
                pictureBox1.Image = bmp;
                Step.Text = step.ToString();
            }
            
        }

        private void DoWork()
        {
            Network network = new Network();
            network.Set();

            for (int i = 0; i < 100000; i++)
            {
            //    Thread.Sleep(2);
                network.SingleStep();
                int percentage = i/1000 ; // Progress
                if(i%20 == 0)
                {
                    BeginInvoke(new Change(Update), percentage, network.GetNodes, network.GetRegion, i);
                }
                    
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            listen.Suspend();
           // Application.Exit();
        }

        private void buttonStart_Click(object sender, EventArgs e)
        {
            buttonStart.Enabled = false;
            listen = new Thread(new ThreadStart(DoWork));
   //         listen.IsBackground = true;
            listen.Start();

        }

        private void resume_Click(object sender, EventArgs e)
        {
             listen.Resume();
        }

        private void Stop_Click(object sender, EventArgs e)
        {
            if(listen.ThreadState == ThreadState.Running)
                listen.Abort();
            buttonStart.Enabled = true;
        }
    }
}
