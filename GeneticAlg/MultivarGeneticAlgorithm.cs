using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GeneticAlg.Extensions;

namespace GeneticAlg
{
    public class MultivarGeneticAlgorithm
    {
        private static Random _rand = new Random();

        #region Problem parameters

        /// <summary>
        /// Left bound of the search interval
        /// </summary>
        private double[] _a;
        /// <summary>
        /// Right bound of the search interval
        /// </summary>
        private double[] _b;

        /// <summary>
        /// Function to minimize
        /// </summary>
        private Func<double[], double> _f;

        /// <summary>
        /// Number of variables in function
        /// </summary>
        private int _varNum;

        /// <summary>
        /// Precision in significant digits
        /// </summary>
        private int _eps;

        #endregion

        #region Algorithm parameters

        /// <summary>
        /// Size of the individuals
        /// </summary>
        private int[] _n;
        /// <summary>
        /// Length of all individual bitwords
        /// </summary>
        private int _m;
        /// <summary>
        /// Size of population
        /// </summary>
        private int _N_p;
        /// <summary>
        /// Size of tournament
        /// </summary>
        private int _N_t;
        /// <summary>
        /// Probability of mutation
        /// </summary>
        private double _p_m;
        /// <summary>
        /// Probability of crossover
        /// </summary>
        private double _p_c;
        /// <summary>
        /// Maximum unchanging populations halt condition
        /// </summary>
        private int _t_max_i;
        /// <summary>
        /// Maximum overall populations halt condition
        /// </summary>
        private int _t_max;

        #endregion

        #region Algorithm variables

        private int t;

        private List<double> _popFitnessStack;

        #endregion


        #region Algorithm functions

        /// <summary>
        /// Return the number of bits for the individual that is required for specified precision
        /// </summary>
        /// <returns></returns>
        private int find_indiv_size(double a_i, double b_i)
        {
            int numOfPoints = (int)Math.Floor((b_i - a_i) * Math.Pow(10, _eps));
            int count = 0;

            while (numOfPoints > 0)
            {
                numOfPoints >>= 1;
                count++;
            }

            return count;
        }

        /// <summary>
        /// Converts a bit word to an integer
        /// </summary>
        /// <param name="v_i"></param>
        /// <returns></returns>
        private double b2d(string v_i, int i)
        {
            int v_i_int = Convert.ToInt32(v_i, 2);
            return _a[i] + (((_b[i] - _a[i]) / (Math.Pow(2, _n[i]) - 1)) * v_i_int);
        }

        /// <summary>
        /// Splits the individual bitword into separate real values
        /// </summary>
        /// <param name="v_i"></param>
        /// <returns></returns>
        private double[] indiv_to_double_array(string v_i)
        {
            var param = new double[_varNum];
            int v_i_idx = 0;

            for (int i = 0; i < _varNum; i++)
            {
                param[i] = b2d(v_i.Substring(v_i_idx, _n[i]), i);
                v_i_idx += _n[i];
            }

            return param;
        }

        /// <summary>
        /// Determines how close to the solution the individual is
        /// </summary>
        /// <param name="v_i"></param>
        /// <returns></returns>
        private double fitness(string v_i)
        {
            return _f(indiv_to_double_array(v_i));
        }

       

        /// <summary>
        /// Determines the fitness of the current population
        /// </summary>
        /// <returns></returns>
        private (double val, double[] x) fitness_pop(List<string> pop)
        {
            var min = pop.Min(fitness);
            return (min, indiv_to_double_array(pop.Find(ind => fitness(ind) <= min)));
        }

        private List<string> tournament(List<string> pop)
        {
            var selectedPop = new List<string>(_N_p);
            var selection = new List<string>(_N_t);

            while (selectedPop.Count < _N_p)
            {
                selection.Clear();

                // pick N_t random elements from current population
                for (int i = 0; i < _N_t; i++)
                {
                    selection.Add(pop[_rand.Next(_N_p)]);
                }

                //Console.WriteLine("Tournament: ");
                //foreach (var sel in selection)
                //{
                //    Console.WriteLine($"\t{sel}");
                //}

                var min = selection.Min(fitness);
                selectedPop.Add(selection.First(v => fitness(v) <= min));
            }

            return selectedPop;
        }

        private List<string> mutate(List<string> pop)
        {
            var mutatedPop = new List<string>(_N_p);

            //Console.WriteLine("Mutation: ");

            foreach (var v in pop)
            {
                if (_rand.NextDouble() <= _p_m)
                {
                    var i_m = _rand.Next(_m);
                    var mutatedInd = v.ToCharArray();

                    mutatedInd[i_m] = (mutatedInd[i_m] == '0') ? '1' : '0';
                    // flip i-th bit
                    mutatedPop.Add(new string(mutatedInd));

                    //Console.WriteLine($"{v} mutated in bit {i_m}");
                }
                else
                {
                    mutatedPop.Add(v);
                }
            }

            return mutatedPop;
        }

        private List<string> crossover(List<string> pop)
        {
            var crossedOverPop = new List<string>(_N_p);

            //Console.WriteLine("Crossover: ");

            while (crossedOverPop.Count < _N_p)
            {
                var v_i = pop[_rand.Next(_N_p)]; // e.g. 01001000
                var v_j = pop[_rand.Next(_N_p)]; // e.g. 11100001

                if (_rand.NextDouble() < _p_c)
                {
                    var i_c = _rand.Next(_m); // e.g. 5

                    var v_i_first_half = v_i.Substring(0, v_i.Length - i_c); // 01001000 & 11100000 = 01000000
                    var v_i_second_half = v_i.Substring(v_i.Length - i_c); // 01001000 & 00011111 = 00001000

                    var v_j_first_half = v_j.Substring(0, v_j.Length - i_c); // 11100001 & 11100000 = 11100000
                    var v_j_second_half = v_j.Substring(v_j.Length - i_c); // 11100001 & 00011111 = 00000001

                    crossedOverPop.Add(string.Concat(v_i_first_half, v_j_second_half)); // 01000000 + 00000001 = 01000001
                    crossedOverPop.Add(string.Concat(v_j_first_half, v_i_second_half)); // 11100000 + 00001000 = 11101000
                }
                else // leave as is
                {
                    crossedOverPop.Add(v_i);
                    crossedOverPop.Add(v_j);
                }
            }

            return crossedOverPop;
        }

        #endregion

        public MultivarGeneticAlgorithm(
            double[] a,
            double[] b,
            Func<double[], double> f,
            int varNum,
            int eps,

            int N_p = 100,
            int N_t = 2,
            double p_m = 0.01,
            double p_c = 0.8,
            int t_max_i = 10,
            int t_max = 1000)
        {
            _a = a;
            _b = b;

            _f = f;
            _varNum = varNum;

            _eps = eps;

            _N_p = N_p;
            
            if (N_t < 2) { throw new ArgumentOutOfRangeException("N_t must be at least 2"); }
            _N_t = N_t; 

            if (p_m < 0 || p_m > 1) { throw new ArgumentOutOfRangeException("p_m must be from 0 to 1"); } 
            _p_m = p_m;

            if (p_c < 0 || p_c > 1) { throw new ArgumentOutOfRangeException("p_c must be from 0 to 1"); }
            _p_c = p_c;

            _t_max_i = t_max_i;
            _t_max = t_max;

            _n = new int[_varNum];
            for (int i = 0; i < _varNum; i++)
            {
                if (_a[i] >= _b[i]) { throw new ArgumentException("a must be smaller than b"); }

                _n[i] = find_indiv_size(_a[i], _b[i]);
                if (_n[i] > 31) { throw new ArgumentException("Precision is too high"); }

                _m += _n[i];
            }
        }


        public double Solve()
        {
            var pop = new List<string>();
            // generate initial population
            for (int i = 0; i < _N_p; i++)
            {
                var indiv = string.Empty;

                for (int j = 0; j < _varNum; j++)
                {
                    var indiv_int = Convert.ToString((uint)_rand.Next((int) Math.Pow(2, _n[j])), 2);
                    //Console.WriteLine(indiv_int);
                    indiv_int = indiv_int.PadLeft(_n[j], '0');

                    indiv += indiv_int;
                }

                pop.Add(indiv);
            }

            int t_i = 0;

            t = 0;
            _popFitnessStack = new List<double>();

            while (t < _t_max)
            {
                var selectedPop = tournament(pop);
                var mutatedPop = mutate(selectedPop);
                var crossedOverPop = crossover(mutatedPop);

                pop = crossedOverPop;

                var fitness = fitness_pop(pop);
                _popFitnessStack.Add(fitness.val);

                Console.WriteLine($"Iteration #{t}: {fitness.val} at x = [{fitness.x.ToArrayString()}]");

                if (t > 1)
                {
                    if (Math.Abs(_popFitnessStack[t] - _popFitnessStack[t - 1]) < Math.Pow(10, -_eps))
                    {
                        t_i++;
                    }
                }

                if (t_i >= _t_max_i)
                {
                    break;
                }

                t++;
            }

            return _popFitnessStack.Last();
        }
    }
}
