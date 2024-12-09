
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using UnityEngine;
using Random = System.Random;

public class HashTable
{
    private const int MaxBucketSize = 512; // Maximum items per bucket
    private const int NumBuckets = 409; // Number of buckets
    private const int MaxIterations = 25; // Maximum retries for cuckoo hashing
    private const int PrimeNumber = 1900813; // A large prime number
    private readonly int[] start;
    private readonly int[] count;
    private readonly ConcurrentDictionary<int, int> uniqueKeys; // Tracks unique keys
    private readonly int[] globalStorage; // Stores values for keys contiguously
    private readonly ConcurrentDictionary<int, Tuple<int, int>> keyIndexCountMap; // Maps key to index and count

    public HashTable(int capacity)
    {
        int totalBuckets = (int)Math.Ceiling((double)capacity / NumBuckets);
        start = new int[totalBuckets];
        count = new int[totalBuckets];
        uniqueKeys = new ConcurrentDictionary<int, int>();
        globalStorage = new int[capacity * 2]; // Storage for values
        keyIndexCountMap = new ConcurrentDictionary<int, Tuple<int, int>>();
    }

    private int Hash1(int key, int bucketCount) => key % bucketCount;

    private int Hash2(int key, int bucketCount, int c0, int c1) =>
        (int)((c0 + (long)c1 * key) % PrimeNumber % bucketCount);

    public void Build(KeyValuePair<int, int>[] items)
    {
        // Phase 1: Distribute items into buckets
        Parallel.For(0, items.Length, i =>
        {
            int bucketIndex = Hash1(items[i].Key, start.Length);
            lock (count) // Atomic increment simulation
            {
                if (count[bucketIndex] < MaxBucketSize)
                {
                    int offset = start[bucketIndex] + count[bucketIndex]++;
                    globalStorage[offset] = items[i].Value;
                    uniqueKeys.AddOrUpdate(items[i].Key, 1, (key, existing) => existing + 1);
                }
                else
                {
                    // Handle rehashing if bucket overflows
                    int c0 = new Random().Next();
                    int c1 = new Random().Next();
                    bucketIndex = Hash2(items[i].Key, start.Length, c0, c1);
                    int offset = start[bucketIndex] + count[bucketIndex]++;
                    globalStorage[offset] = items[i].Value;
                    uniqueKeys.AddOrUpdate(items[i].Key, 1, (key, existing) => existing + 1);
                }
            }
        });

        // Phase 2: Build cuckoo hash tables for each bucket
        Parallel.For(0, start.Length, bucketIndex =>
        {
            var localItems = uniqueKeys.Where(kv => Hash1(kv.Key, start.Length) == bucketIndex).ToArray();
            BuildCuckooHashTable(localItems, bucketIndex);
        });
    }

    private void BuildCuckooHashTable(KeyValuePair<int, int>[] items, int bucketIndex)
    {
        // Initialize sub-tables
        var T1 = new int[192];
        var T2 = new int[192];
        var T3 = new int[192];

        int c0 = new Random().Next();
        int c1 = new Random().Next();

        for (int iteration = 0; iteration < MaxIterations; iteration++)
        {
            bool allPlaced = true;

            foreach (var item in items)
            {
                int hash1 = Hash1(item.Key, T1.Length);
                if (T1[hash1] == 0)
                {
                    T1[hash1] = item.Key;
                    keyIndexCountMap[item.Key] = new Tuple<int, int>(hash1, item.Value);
                    continue;
                }

                int hash2 = Hash2(item.Key, T2.Length, c0, c1);
                if (T2[hash2] == 0)
                {
                    T2[hash2] = item.Key;
                    keyIndexCountMap[item.Key] = new Tuple<int, int>(hash2, item.Value);
                    continue;
                }

                int hash3 = Hash2(item.Key, T3.Length, c0 + 1, c1 + 1);
                if (T3[hash3] == 0)
                {
                    T3[hash3] = item.Key;
                    keyIndexCountMap[item.Key] = new Tuple<int, int>(hash3, item.Value);
                    continue;
                }

                // If not placed, choose another random hash and retry
                allPlaced = false;
                c0 = new Random().Next();
                c1 = new Random().Next();
                break;
            }

            if (allPlaced)
                break;
        }

        // Rearrange data for contiguous storage
        int globalOffset = start[bucketIndex];
        foreach (var kv in keyIndexCountMap)
        {
            int index = kv.Value.Item1;
            int count = kv.Value.Item2;
            Array.Copy(globalStorage, index, globalStorage, globalOffset, count);
            globalOffset += count;
        }
    }

    public int[] Retrieve(int key)
    {
        if (keyIndexCountMap.TryGetValue(key, out var data))
        {
            int index = data.Item1;
            int count = data.Item2;
            var values = new int[count];
            Array.Copy(globalStorage, index, values, 0, count);
            return values;
        }

        throw new Exception("Key not found.");
    }
}

