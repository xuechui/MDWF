using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MDtestwinform
{
    public class TBuf   //time-dependent property
    {
        private List<VecR> orgR = new List<VecR>();
        private List<VecR> rTrue = new List<VecR>();
        List<double> rrDiffuse = new List<double>();
        private int count = 0;

        public int GetCount { get { return count; } }
        public void SetCount(int a) { count = a; }
        public List<VecR> GetorgR { get { return orgR; } }
        public List<VecR> GetrTrue { get { return rTrue; } }
        public List<double> GetrrDiffuse { get { return rrDiffuse; } }


        public void AddorigR(VecR r)
        {
            orgR.Add(r);
        }
        public void AddrTrue(VecR r)
        {
            rTrue.Add(r);
        }
        public void InitrrDiffuse(double r)
        {
            rrDiffuse.Add(r);
        }
        public void AddrrDiffuse(int i,double r)
        {
            rrDiffuse[i] += r;
        }
        public void AddCount()
        {
            count++;
        }

    }

    public struct VecR
    {
        public double x, y, z;
        public VecR(double x, double y, double z)
        {
            this.x = x; this.y = y; this.z = z;
        }
        public VecR Get { get { return this; } }

        public void VSub(VecR a, VecR b)
        {
            x = a.x + b.x; y = a.y + b.y; z = a.z + b.z;
        }
        public void VDiv(VecR a, VecR b)
        {
            x = a.x / b.x; y = a.y / b.y; z = a.z / b.z;
        }
        public void VMul(VecR a, VecR b)
        {
            x = a.x * b.x; y = a.y * b.y; z = a.z * b.z;
        }
        public void VAdd(VecR a, VecR b)
        {
            x = a.x + b.x; y = a.y + b.y; z = a.z + b.z;
        }

        public void VZero()
        { x = 0; y = 0; z = 0; }
        public void VWrap(VecR region)
        {
            if (x >= 0.5 * region.x) { x -= region.x; }
            if (x < -0.5 * region.x) { x += region.x; }
            if (y >= 0.5 * region.y) { y -= region.y; }
            if (y < -0.5 * region.y) { y += region.y; }
            if (z >= 0.5 * region.z) { z -= region.z; }
            if (z < -0.5 * region.z) { z += region.z; }
        }
        public void VVAdd(VecR r)
        {
            x += r.x; y += r.y; z += r.z;
        }
        public void Nint()  //Return Round-up Int
        {
            x = (x < 0) ? (-(int)(0.5 - x)) : ((int)(0.5 + x));
            y = (y < 0) ? (-(int)(0.5 - y)) : ((int)(0.5 + y));
            z = (z < 0) ? (-(int)(0.5 - z)) : ((int)(0.5 + z));
        }
        public double VLenSq()
        {
            return x * x + y * y + z * z;
        }
    }

    public struct VecI
    {
        public int x, y, z;
        public VecI(int x, int y, int z)
        {
            this.x = x; this.y = y; this.z = z;
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

        public Node(VecR r, VecR rv, VecR ra)
        {
            this.r = r;
            this.rv = rv;
            this.ra = ra;
        }
        public VecR Getr { get { return r; } }
        public VecR Getv { get { return rv; } }
        public void Setr(VecR r) { this.r = r; }
        public void VZeroA()     //Set Acceleration to zero
        { ra.x = 0; ra.y = 0; ra.z = 0; }

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


    public struct Prop
    {
        double val, sum, sum2;
        public Prop(double val, double sum, double sum2)
        {
            this.val = val;
            this.sum = sum;
            this.sum2 = sum2;
        }
        public void SetVal(double value)
        {
            val = value;
        }
        public double GetVal { get { return val; } }

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
            if (sum2 <= 0)
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
