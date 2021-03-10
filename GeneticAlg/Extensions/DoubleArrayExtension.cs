using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeneticAlg.Extensions
{
    public static class DoubleArrayExtension
    {
        public static string ToArrayString(this double[] doubleArray)
        {
            var res = string.Empty;
            foreach (var d in doubleArray)
            {
                res += $"{d}, ";
            }

            return res.Substring(0, res.Length - 2);
        }
    }
}
