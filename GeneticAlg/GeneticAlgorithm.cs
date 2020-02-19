using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeneticAlg
{
    public class GeneticAlgorithm
    {
        private static Random _rand = new Random();

        #region Problem parameters

        /// <summary>
        /// Left bound of the search interval
        /// </summary>
        private double _a;
        /// <summary>
        /// Right bound of the search interval
        /// </summary>
        private double _b;

        /// <summary>
        /// Function to minimize
        /// </summary>
        private Func<double, double> _f;

        /// <summary>
        /// Precision in significant digits
        /// </summary>
        private int _eps;

        #endregion

        #region Algorithm parameters

        /// <summary>
        /// Size of the individual
        /// </summary>
        private int _n;
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
        private int find_indiv_size()
        {
            int numOfPoints = (int)Math.Floor((_b - _a) * Math.Pow(10, _eps));
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
        private double b2d(int v_i)
        {
            return _a + (((_b - _a) / (Math.Pow(2, _n) - 1)) * v_i);
        }

        /// <summary>
        /// Determines how close to the solution the individual is
        /// </summary>
        /// <param name="v_i"></param>
        /// <returns></returns>
        private double fitness(int v_i)
        {
            return _f(b2d(v_i));
        }

        /// <summary>
        /// Determines the fitness of the current population
        /// </summary>
        /// <returns></returns>
        private (double val, double x) fitness_pop(List<int> pop)
        {
            var min = pop.Min(fitness);
            return (min, b2d(pop.Find(ind => fitness(ind) <= min)));
        }

        private List<int> tournament(List<int> pop)
        {
            var selectedPop = new List<int>(_N_p);
            var selection = new List<int>(_N_t);

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
                //    Console.WriteLine($"\t{Convert.ToString(sel, 2)}");
                //}

                var min = selection.Min(fitness);
                selectedPop.Add(selection.First(v => fitness(v) <= min));
            }

            return selectedPop;
        }

        private List<int> mutate(List<int> pop)
        {
            var mutatedPop = new List<int>(_N_p);

            //Console.WriteLine("Mutation: ");

            foreach (var v in pop)
            {
                if (_rand.NextDouble() <= _p_m)
                {
                    var i_m = _rand.Next(_n);
                    var mask = 1 << i_m;

                    // flip i-th bit
                    mutatedPop.Add(v ^ mask);

                    //Console.WriteLine($"{Convert.ToString(v, 2)} mutated in bit {i_m}");
                }
                else
                {
                    mutatedPop.Add(v);
                }
            }

            return mutatedPop;
        }

        private List<int> crossover(List<int> pop)
        {
            var mutatedPop = new List<int>(_N_p);

            //Console.WriteLine("Crossover: ");

            while (mutatedPop.Count < _N_p)
            {
                int v_i = pop[_rand.Next(_N_p)]; // e.g. 01001000
                int v_j = pop[_rand.Next(_N_p)]; // e.g. 11100001

                if (_rand.NextDouble() < _p_c)
                {
                    var i_c = _rand.Next(_n); // e.g. 5

                    var v_i_first_half = v_i & (int.MaxValue << i_c); // 01001000 & 11100000 = 01000000
                    var v_i_second_half = v_i & ((int)Math.Pow(2, i_c) - 1); // 01001000 & 00011111 = 00001000

                    var v_j_first_half = v_j & (int.MaxValue << i_c); // 11100001 & 11100000 = 11100000
                    var v_j_second_half = v_j & ((int)Math.Pow(2, i_c) - 1); // 11100001 & 00011111 = 00000001

                    mutatedPop.Add(v_i_first_half + v_j_second_half); // 01000000 + 00000001 = 01000001
                    mutatedPop.Add(v_j_first_half + v_i_second_half); // 11100000 + 00001000 = 11101000

                    //Console.WriteLine($"\nv_i: {Convert.ToString(v_i, 2)}");
                    //Console.WriteLine($"v_j: {Convert.ToString(v_j, 2)}");
                    
                    //Console.WriteLine($"i_c: {i_c}");

                    //Console.WriteLine($"v_i_first_half: {Convert.ToString(v_i_first_half, 2)}");
                    //Console.WriteLine($"v_i_second_half: {Convert.ToString(v_i_second_half, 2)}");
                    //Console.WriteLine($"v_j_first_half: {Convert.ToString(v_j_first_half, 2)}");
                    //Console.WriteLine($"v_j_second_half: {Convert.ToString(v_j_second_half, 2)}");


                }
                else // leave as is
                {
                    mutatedPop.Add(v_i);
                    mutatedPop.Add(v_j);
                }
            }

            return mutatedPop;
        }

        #endregion

        public GeneticAlgorithm(
            double a,
            double b,
            Func<double, double> f,
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

            if (_a >= _b) { throw new ArgumentException("a must be smaller than b"); }

            _f = f;

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

            _n = find_indiv_size();
            if (_n > 31) { throw new ArgumentException("Precision is too high"); }
        }


        public double Solve()
        {
            var pop = new List<int>();
            // generate initial population
            for (int i = 0; i < _N_p; i++)
            {
                pop.Add(_rand.Next((int)Math.Pow(2, _n)));
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

                //Console.WriteLine($"Iteration #{t}: {fitness.val} at x = {fitness.x}");

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
