using System;

namespace GeneticAlg
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Greets!");

            var solver = new GeneticAlgorithm(
                -1,
                2,
                x => (x * Math.Sin(10 * Math.PI * x)) + 1,
                6,
                100,
                2);

            Console.WriteLine($"\nSolution: {solver.Solve()}");

            Console.ReadKey();
            Console.WriteLine();

            var popNums = new int[6] { 10, 50, 100, 200, 500, 1000 };

            foreach (var num in popNums)
            {
                solver = new GeneticAlgorithm(
                    -1,
                    2,
                    x => (x * Math.Sin(10 * Math.PI * x)) + 1,
                    6,
                    num,
                    2);

                Console.WriteLine($"Solution for population {num}: {solver.Solve()}");
            }

            Console.ReadKey();
            Console.WriteLine();

            var P_cs = new double[5] { 1, 0.9, 0.8, 0.5, 0.1 };

            foreach (var num in P_cs)
            {
                solver = new GeneticAlgorithm(
                    -1,
                    2,
                    x => (x * Math.Sin(10 * Math.PI * x)) + 1,
                    6,
                    100,
                    2,
                    p_c: num);

                Console.WriteLine($"Solution for crossover prob {num}: {solver.Solve()}");
            }

            Console.ReadKey();
            Console.WriteLine();

            var N_ts = new int[5] { 2, 3, 4, 5, 6 };

            foreach (var num in N_ts)
            {
                solver = new GeneticAlgorithm(
                    -1,
                    2,
                    x => (x * Math.Sin(10 * Math.PI * x)) + 1,
                    6,
                    100,
                    num);

                Console.WriteLine($"Solution for N_t {num}: {solver.Solve()}");
            }

            Console.ReadKey();
        }
    }
}
