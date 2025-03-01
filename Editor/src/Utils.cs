
using System;
using System.Threading;

namespace TerraMesh.Utils
{
    /// <summary>
    /// A thread-safe random number generator.
    /// </summary>
    /// <remarks>
    /// This class is based on the implementation by Andrew Lock:
    /// https://andrewlock.net/building-a-thread-safe-random-implementation-for-dotnet-framework/
    /// </remarks>
    internal class ThreadSafeRandom : Random
    {
        [ThreadStatic]
        private static Random? _local;

        private static readonly Random _global = new Random(); // used for thread-safe seed randomization

        // Constructor to allow seeding from external sources
        public ThreadSafeRandom() : this(GenerateSeed())
        {
        }

        public ThreadSafeRandom(int seed) : base(seed)
        {
            // We don't actually use the seed passed to the base constructor directly,
            // but we call the base constructor to fulfill the inheritance requirement.
            // The thread-local instance will handle the actual random number generation.
            _local = new Random(seed);
        }

        private static int GenerateSeed()
        {
            lock (_global)
            {
                return _global.Next();
            }
        }

        internal static Random Instance
        {
            get
            {
                if (_local is null)
                {
                    _local = new Random(GenerateSeed());
                }
                return _local;
            }
        }

        public override int Next() => Instance.Next();

        public override int Next(int maxValue) => Instance.Next(maxValue);

        public override int Next(int minValue, int maxValue) => Instance.Next(minValue, maxValue);

        public override void NextBytes(byte[] buffer) => Instance.NextBytes(buffer);
        
        public override double NextDouble() => Instance.NextDouble();
        
    }
}