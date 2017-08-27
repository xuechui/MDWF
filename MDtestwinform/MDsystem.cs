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
        public double GetDen { get { return density; } }
    }


    class Network
    {
        private int nMol;
        private int stepCount = 0;
        private double timeNow = 0;
        private double uSum = 0; private double virSum = 0;
        private VecR vSum;
        private double vvSum;
        private double deltaT = 0.006;
        private int NDIM = 3;  //system dimension
        private double temperature = 1.0;
        private double avgChainLen;

        private List<Node> nodes;
        //Parameters:
        private double rCut, density;
        private VecR region;
        private VecI initUcell;
        private double velMag;   //To scale velcity
        private Prop kinEnergy, pressure, totEnergy;

        /// <summary>
        /// Parameters to calculate diffusion
        /// </summary>
        private int countDiffuseAv, limitDiffusivaAv, nBuffDiffuse,
            nValDiffuse, stepDiffuse;
        List<TBuf> tBuf;
        private List<double> rrDiffuseAv = new List<double>();
        
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

        public void TransferedParams(Param para)
        {
            density = para.GetDen;
        }

        public void SetParams()
        {
            //Cut-off radius
            rCut = Math.Pow(2.0, 1.0 / 6.0);
    //        density = 0.85;
            initUcell = new VecI(3, 3, 3);  //Firstly manually assign the cell
            nMol = initUcell.x * initUcell.y * initUcell.z;
            velMag = Math.Sqrt(NDIM * (1.0 - 1.0 / nMol) * temperature);

            double fac = 1.0 / Math.Pow(density, 1.0 / 3.0);
            region = new VecR(fac * initUcell.x, fac * initUcell.y, fac * initUcell.z);

            //Set Diffusion Params
            nValDiffuse = 10;
            nBuffDiffuse = 1;
            limitDiffusivaAv = 100;
            stepDiffuse = 10;

        }
        public void SetupJob()
        {
            InitCoords(initUcell, density);
            InitVels();
            InitAccels();
            InitDiffusion();
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
            if(stepCount % 1000 == 0)
                Console.WriteLine(avgChainLen + "and" + pressure.GetVal);
            if (stepCount % stepDiffuse == 0)
            {
                EvalDiffusion();
                
            }
               
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
                     //       fcVal = -90.6 / (1.0 - rr / 2.25);
                            fcVal = -20 / (1.0 - rr / 2.25);
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

            avgChainLen = 0;
            double chainLen = 0;
            int chainNum = 0;
            VecR dr;

            foreach (Node node in nodes)
            {
                vSum.VVAdd(node.Getv);
                vv = node.Getv.x * node.Getv.x + node.Getv.y * node.Getv.y + node.Getv.z * node.Getv.z;
                vvSum += vv;
            }
            kinEnergy.SetVal(0.5 * vvSum / nMol);
            totEnergy.SetVal(kinEnergy.GetVal + uSum / nMol);
            pressure.SetVal(density * (vvSum + virSum)/(nMol * NDIM)  );


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
                        chainLen += dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                        chainNum ++;
                    }
                }
            }
            avgChainLen = chainLen / chainNum;
        }

        void EvalDiffusion()
        {
            VecR dr = new VecR();
            int  ni;
            for (int nb = 0; nb < nBuffDiffuse; nb++)
            {
                if (tBuf[nb].GetCount == 0)
                {
                    foreach (Node node in nodes)
                    {
                        tBuf[nb].AddorigR(node.Getr);
                        tBuf[nb].AddrTrue(node.Getr);
                    }
                }
                if (tBuf[nb].GetCount >= 0)
                {
                    ni = tBuf[nb].GetCount;
                    tBuf[nb].InitrrDiffuse(0.0);
                    for(int n = 0; n < nodes.Count; n ++)
                    {
                        dr.VSub(tBuf[nb].GetrTrue[n].Get, nodes[n].Getr);
                        dr.VDiv(dr.Get, region.Get);
                        dr.Nint();
                        dr.VMul(dr.Get, region.Get);
                        tBuf[nb].GetrTrue[n].Get.VAdd(nodes[n].Getr, dr.Get);
                        dr.VSub(tBuf[nb].GetrTrue[n].Get, tBuf[nb].GetorgR[n].Get);
                        tBuf[nb].AddrrDiffuse(ni, dr.VLenSq());
                    }
                }
                tBuf[nb].AddCount();
            }
            AccumDiffusion();
        }

        public void AccumDiffusion()
        {
            double fac;
            for(int nb = 0; nb < nBuffDiffuse; nb ++)
            {
                if(tBuf[nb].GetCount == nValDiffuse)
                {
                    for(int j = 0; j < nValDiffuse; j++)
                        rrDiffuseAv[j] += tBuf[nb].GetrrDiffuse[j];
                    tBuf[nb].SetCount(0);
                    ++countDiffuseAv;
                    if(countDiffuseAv == limitDiffusivaAv)
                    {
                        fac = 1.0 / (NDIM * 2 * nMol * stepDiffuse *
                            deltaT * limitDiffusivaAv);
                        for (int j = 1; j < nValDiffuse; j++)
                            rrDiffuseAv[j] *= fac / j;
                        Console.WriteLine("diffusion:" + rrDiffuseAv[8]);
                        ZeroDiffusion();
                    }
                }
            }
        }

        public void InitDiffusion()
        {
            tBuf = new List<TBuf>();
            for (int nb = 0; nb < nBuffDiffuse; nb++)
            {
                tBuf.Add(new TBuf());
                tBuf[nb].SetCount(-nb * nValDiffuse / nBuffDiffuse);
            }
            countDiffuseAv = 0;
            for (int j = 0; j < nValDiffuse; j++)
                rrDiffuseAv.Add(0);
        }
        public void ZeroDiffusion()
        {
            countDiffuseAv = 0;
            for (int j = 0; j < nValDiffuse; j++)
                rrDiffuseAv[j] = 0;
        }

    }


    

    
}
