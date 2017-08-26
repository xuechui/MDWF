using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MDtestwinform
{
    public struct Param
    {
        public double density;
        public Param(double den)
        {
            density = den;
        }
    }


    class MDsystem
    {
        private Network network;
    }
    class Network
    {
        private int nMol;
        private int stepCount = 0;
        private double timeNow = 0;
        private double uSum = 0; private double virSum = 0;
        private VecR vSum;
        private double vvSum;
        private double deltaT = 0.01;
        private int NDIM = 3;  //system dimension
        private double temperature = 10;

        private List<Node> nodes;
        //Parameters:
        private double rCut, density;
        private VecR region;
        private VecI initUcell;
        private double velMag;   //To scale velcity
        private Prop kinEnergy, pressure, totEnergy;

        public List<Node> GetNodes { get { return nodes; } }

        public VecR GetRegion { get { return region; } }

        /// <summary>
        /// Set Parameters.
        /// </summary>
        public void Set()
        {
            nodes = new List<Node>();
            SetParams();
            SetupJob();
        }

        public void SetParams()
        {
            //Cut-off radius
            rCut = Math.Pow(2.0, 1.0 / 6.0);
            density = 0.8;
            initUcell = new VecI(3, 3, 3);  //Firstly manually assign the cell
            nMol = initUcell.x * initUcell.y * initUcell.z;
            velMag = Math.Sqrt(NDIM * (1.0 - 1.0 / nMol) * temperature);

            double fac = 1.0 / Math.Pow(density, 1.0 / 3.0);
            region = new VecR(fac * initUcell.x, fac * initUcell.y, fac * initUcell.z);
        }
        public void SetupJob()
        {
            InitCoords(initUcell, density);
            InitVels();
            InitAccels();
            Console.WriteLine("Setup job...deltaT = " + deltaT);
        }
        //Initialize Coordinates
        public void InitCoords(VecI initUcell, double density)
        {
            int nx, ny, nz;
            VecR c, gap;
            int n = 0;

            Console.WriteLine("density = " + density);
            gap = new VecR(region.x / initUcell.x, region.y / initUcell.y, region.z / initUcell.z);

            for (nz = 0; nz < initUcell.z; nz++)
            {
                for (ny = 0; ny < initUcell.y; ny++)
                {
                    for (nx = 0; nx < initUcell.x; nx++)
                    {
                        c.x = nx + 0.5; c.y = ny + 0.5; c.z = nz + 1.5;
                        c.x *= gap.x; c.y *= gap.y; c.z *= gap.z;
                        c.x += -0.5 * region.x; c.y += -0.5 * region.y; c.z += -0.5 * region.z;
                        Node b = new Node(c, c, c);
                        nodes.Add(b);
                        //Introduce connection

                        if (nx > 0.5)
                        {
                            nodes[n - 1].AddCon(n);
                            nodes[n].AddCon(n - 1);
                        }
                        n++;
                    }
                }
            }
        }

        public void InitVels()
        {
            vSum.x = 0; vSum.y = 0; vSum.z = 0;
            foreach (Node node in nodes)
            {
                node.VRand();
                node.VScale(velMag);
                vSum.x += node.Getv.x;
                vSum.y += node.Getv.y;
                vSum.z += node.Getv.z;
            }
            foreach (Node node in nodes)
            {
                node.VVSAddV(-1.0 / nMol, vSum);
            }
        }
        public void InitAccels()
        {
            foreach (Node node in nodes)
            {
                node.VZeroA();
            }
        }

        /// <summary>
        /// Define MD procedure.
        /// </summary>

        public void SingleStep()
        {
            ++stepCount;
            timeNow = stepCount * deltaT;
            LeapfrogStep(1);
            ApplyBoundaryCond();
            ComputeForces();
            ComputeFENE();
            AddFriction(0.5);
            AddWhiteNoise();
            LeapfrogStep(2);
            EvalProps();
            Console.WriteLine(pressure.GetVal);
        }

        public void ComputeForces()
        {
            VecR dr;
            double fcVal, rr, rrCut, rri, rri3;
            rrCut = rCut * rCut;
            foreach (Node node in nodes)
            {
                node.VZeroA();
            }
            uSum = 0;
            virSum = 0;
            for (int j1 = 0; j1 < nodes.Count - 1; j1++)
            {
                for (int j2 = j1 + 1; j2 < nodes.Count; j2++)
                {
                    dr.x = nodes[j1].Getr.x - nodes[j2].Getr.x;
                    dr.y = nodes[j1].Getr.y - nodes[j2].Getr.y;
                    dr.z = nodes[j1].Getr.z - nodes[j2].Getr.z;

                    dr.VWrap(region);
                    rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                    if (rr < rrCut)  //within rc
                    {
                        rri = 1.0 / rr;
                        rri3 = rri * rri * rri;
                        fcVal = 48.0 * rri3 * (rri3 - 0.5) * rri;
                        nodes[j1].VVSAddR(fcVal, dr);
                        nodes[j2].VVSAddR(-fcVal, dr);
                        uSum = 4.0 * rri3 * (rri3 - 1.0) + 1.0;
                        virSum += fcVal * rr;
                    }
                }
            }
        }

        public void ComputeFENE()
        {
            VecR dr;
            double fcVal, rr;

            //FENE
            for (int j1 = 0; j1 < nodes.Count; j1++)
            {
                foreach (int j2 in nodes[j1].GetCon)
                {
                    if (j1 < j2)
                    {
                        dr.x = nodes[j1].Getr.x - nodes[j2].Getr.x;
                        dr.y = nodes[j1].Getr.y - nodes[j2].Getr.y;
                        dr.z = nodes[j1].Getr.z - nodes[j2].Getr.z;

                        dr.VWrap(region);
                        rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                        if (rr < 2.25)  //within rc
                        {
                            fcVal = -90.6 / (1.0 - rr / 2.25);
                            nodes[j1].VVSAddR(fcVal, dr);
                            nodes[j2].VVSAddR(-fcVal, dr);
                        }
                    }
                }
            }
        }
        public void AddFriction(double per)
        {
            foreach (Node node in nodes)
            {
                node.VVSAddR(-1.0 * per, node.Getv);
            }
        }
        public double GaussianRand()
        {
            double V1, V2, S;
            double X = 0;
            Random d = new Random((int)DateTime.Now.Ticks);

            do
            {
                double U1 = d.NextDouble();
                double U2 = d.NextDouble();

                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;

                S = V1 * V1 + V2 * V2;
            } while (S >= 1 || S == 0);
            X = V1 * Math.Sqrt(-2 * Math.Log(S) / S);

            return X;
        }

        public void AddWhiteNoise()
        {
            foreach(Node node in nodes)
            {
                node.Noise(GaussianRand() * temperature * 6.0);
            }
            
        }


            public void LeapfrogStep(int part)
        {
            if(part == 1)
            {
                foreach(Node node in nodes)
                {
                    node.LeapFrog1(deltaT);
                }
            }
            else
            {
                foreach(Node node in nodes)
                {
                    node.LeapFrog2(deltaT);
                }
            }
        }
        public void ApplyBoundaryCond()
        {
            foreach(Node node in nodes)
            {
                node.VWrap(region);
            }
        }
        public void AccumProps(int icode)
        {

        }
        public void EvalProps()
        {
            double vv;
            vSum.VZero();
            vvSum = 0;
            foreach (Node node in nodes)
            {
                vSum.VVAdd(node.Getv);
                vv = node.Getv.x * node.Getv.x + node.Getv.y * node.Getv.y + node.Getv.z * node.Getv.z;
                vvSum += vv;
            }
            kinEnergy.SetVal(0.5 * vvSum / nMol);
            totEnergy.SetVal(kinEnergy.GetVal + uSum / nMol);
            pressure.SetVal(density * (vvSum + virSum)/(nMol * NDIM)  );
        }

    }


    class Node
    {
        private VecR r, rv, ra;   //Coordiate, Velocity, Acceleration
        private List<int> ConNode = new List<int>();

        public List<int> GetCon { get { return ConNode; } }
        public void AddCon(int con)
        {
            ConNode.Add(con);
        }

        public Node (VecR r, VecR rv, VecR ra)
        {
            this.r = r;
            this.rv = rv;
            this.ra = ra;
        }
        public VecR Getr{ get { return r; } }
        public VecR Getv{ get { return rv; } }
        public void Setr(VecR r)  { this.r = r; }
        public void VZeroA()     //Set Acceleration to zero
        { ra.x = 0; ra.y = 0; ra.z = 0;   }

        public void Noise(double noise)
        {
            ra.x += noise;
            ra.y += noise;
            ra.z += noise;
        }

        public void VVSAddR(double fac, VecR dr)
        {
            ra.x += fac * dr.x;
            ra.y += fac * dr.y;
            ra.z += fac * dr.z;
        }
        public void VVSAddV(double fac, VecR dr)
        {
            rv.x += fac * dr.x;
            rv.y += fac * dr.y;
            rv.z += fac * dr.z;
        }
        public void VWrap(VecR region)
        {
            if (r.x >= 0.5 * region.x) { r.x -= region.x; }
            if (r.x < -0.5 * region.x) { r.x += region.x; }
            if (r.y >= 0.5 * region.y) { r.y -= region.y; }
            if (r.y < -0.5 * region.y) { r.y += region.y; }
            if (r.z >= 0.5 * region.z) { r.z -= region.z; }
            if (r.z < -0.5 * region.z) { r.z += region.z; }
        }
        public void VRand()
        {
            Random d = new Random((int)DateTime.Now.Ticks);
            rv.x = d.NextDouble();
    //        Console.WriteLine(rv.x);
            rv.y = d.NextDouble();
            rv.z = d.NextDouble();
        }
        public void VScale(double velMag)
        {
            rv.x *= velMag; rv.y *= velMag; rv.z *= velMag;
        }


        public void LeapFrog1(double deltaT)
        {
            rv.x += 0.5 * deltaT * ra.x; rv.y += 0.5 * deltaT * ra.y; rv.z += 0.5 * deltaT * ra.z;
            r.x += deltaT * rv.x; r.y += deltaT * rv.y; r.z += deltaT * rv.z;
        }
        public void LeapFrog2(double deltaT)
        {
            rv.x += 0.5 * deltaT * ra.x; rv.y += 0.5 * deltaT * ra.y; rv.z += 0.5 * deltaT * ra.z;
        }
    }

    public struct VecR
    {
        public double x, y, z;
        public VecR(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public void VZero()
        { x = 0; y = 0; z = 0; }
        public void VWrap(VecR region)
        {
            if(x >= 0.5 * region.x){ x -= region.x; }
            if(x < -0.5 * region.x){ x += region.x; }
            if (y >= 0.5 * region.y) { y -= region.y; }
            if (y < -0.5 * region.y) { y += region.y; }
            if (z >= 0.5 * region.z) { z -= region.z; }
            if (z < -0.5 * region.z) { z += region.z; }
        }
        public void VVAdd(VecR r)
        {
            x += r.x;
            y += r.y;
            z += r.z;
        }
        

    }

    public struct VecI
    {
        public int x, y, z;
        public VecI(int x, int y, int z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }
    public struct Prop
    {
        double val, sum, sum2;
        public Prop(double val,double sum,double sum2)
        {
            this.val = val;
            this.sum = sum;
            this.sum2 = sum2;
        }
        public void SetVal(double value)
        {
            val = value;
        }
        public double GetVal{ get { return val; } }

        public void PropZero()
        {
            sum = 0;
            sum2 = 0;
        }
        public void PropAccum()
        {
            sum += val;
            sum2 += val;
        }
        public void PropAvg(int n)
        {
            sum /= n;

            sum2 = sum2 / n - Math.Sqrt(sum);
            if(sum2 <= 0)
            {
                sum2 = 0;
            }
            else
            {
                sum2 = Math.Sqrt(sum2);
            }
        }
    }

}
