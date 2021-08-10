using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Windows.Forms;

namespace ImageQuantization
{
    public class edgee {
        public int src, dest;
        public float weight;
        public edgee(int src,int dest,float weight) {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }
    }
    struct Edge : IComparable<Edge>
    {
        // Each Edge Consists Of A Source, Destination, And Weight
        public RGBPixel src, dest;
        public float weight;

        // Constructor
        public Edge(RGBPixel src, RGBPixel dest, float weight)
        {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }

        // In The MST, Edges Have To Be Sorted Ascendingly
        // According To Their Weights, So This
        // Function Is Used To Base The Sorting Operation
        // On Weights
        public int CompareTo(Edge compareEdge)
        {
            if (this.weight > compareEdge.weight)
            {
                return 1;
            }

            else if (this.weight < compareEdge.weight)
            {
                return -1;
            }

            else
                return 0;
        }
    }

    // This Class Represents Each Node In The Priority Queue
    // Each Node Contains A Vertex (Color) And A Pointer
    // To That Node's Parent Vertex
    class PQNode : FastPriorityQueueNode
    {
        public RGBPixel currentVertex;
        public RGBPixel? Parent;
        public PQNode(RGBPixel vertex, RGBPixel? parent)
        {
            currentVertex = vertex;
            Parent = parent;
        }


    }


    class Graph
    {
        // Member Variables
        int vertex;
        public List<Edge> MST;
        public float MSTSum;
        public static Dictionary<int, List<int>> neighbours;
        public static List<edgee> nMST;


        // Constructor
        public Graph(int v)
        {
            this.vertex = v; ;
            MST = new List<Edge>();
            neighbours = new Dictionary<int, List<int>>();
            MSTSum = 0;
            nMST = new List<edgee>();
        }


        public List<edgee> newMST(List<RGBPixel> DOC) {
            List<float> costs = Enumerable.Repeat(float.MaxValue, DOC.Count).ToList();
            List<int> parents = Enumerable.Repeat(-1, DOC.Count).ToList();
            int nextNode = 0;
            bool[] MSTvisted = new bool[DOC.Count];
            RGBPixel firstPixel, secondPixel;
            float cost = 0;
            
            for (int i = 0; i < DOC.Count-1; i++) {
                int index = nextNode;
                firstPixel = DOC[index];
                MSTvisted[index] = true;
                float mininode = float.MaxValue;
                for (int j = 0; j < DOC.Count; j++) {
                    if (!MSTvisted[j])
                    {
                        secondPixel = DOC[j];
                        cost =  (float)Math.Sqrt((firstPixel.red - secondPixel.red) * (firstPixel.red- secondPixel.red) 
                               +(firstPixel.green - secondPixel.green) * (firstPixel.green - secondPixel.green) 
                               +(firstPixel.blue - secondPixel.blue) * (firstPixel.blue - secondPixel.blue));
                        if (cost < costs[j]) {
                            parents[j] = index;
                            costs[j] = cost;
                        }

                        if (costs[j] < mininode) {
                            mininode = costs[j];
                            nextNode = j;
                        }

                    }
                }
               
            }
            MSTSum = 0;
            for (int i = 1; i < parents.Count; i++){
                nMST.Add(new edgee(parents[i],i,costs[i]));
                MSTSum += costs[i];
            }
            //foreach (var i in nMST) {
            //    MessageBox.Show("" + i.weight);
            //}
            //MessageBox.Show(""+MSTSum);
            return nMST;
        }
        #region Spanning Tree Methods

        // Used To Calculate The Weight Difference Between Two Nodes
        private static float WeightDifference(PQNode N1, PQNode N2)
        {
            byte red1, green1, blue1;
            byte red2, green2, blue2;

            // Assigning The RGB Values Of Each Color
            red1 = N1.currentVertex.red;
            red2 = N2.currentVertex.red;
            green1 = N1.currentVertex.green;
            green2 = N2.currentVertex.green;
            blue1 = N1.currentVertex.blue;
            blue2 = N2.currentVertex.blue;

            // Calculating The Difference Between Them
            return (float)Math.Sqrt((red2 - red1) * (red2 - red1) + (green2 - green1) * (green2 - green1) + (blue2 - blue1) * (blue2 - blue1));
        }

        // The Function That Calculates The Minimum Spanning Tree For The Graph
        // Its Only Parameter Is The List Of Distinct Colors
        public void MinSpanningTree(List<RGBPixel> LoD)
        {
            // Using A Priority Queue So That The Edges Are
            // Sorted According To Their Weights In Ascending Order
            // The Priority Queue Contains The Number Of Vertices (LoD.Count)
            FastPriorityQueue<PQNode> PQ = new FastPriorityQueue<PQNode>(LoD.Count);

            // An Array That Contains Each Node (color)
            // And A Pointer To That Node's Parent
            PQNode[] PQNode = new PQNode[LoD.Count];

            // Initialising The First Node In The Array
            // Its Parent Is Null
            PQNode[0] = new PQNode(LoD[0], null);
            // Inserting That First Node Into The Priority Queue
            // With A Value of 0
            PQ.Enqueue(PQNode[0], 0);

            float weight;

            // For Each Vertex (Minus The Already Added One), 
            // Insert That Vertex Into The Array
            // With Parent = null, Then Insert It Into The PQ
            // With Values = Infinity
            for (int i = 1; i < LoD.Count; i++)
            {
                PQNode[i] = new PQNode(LoD[i], null);
                PQ.Enqueue(PQNode[i], int.MaxValue);
            }

            // Looping Until The Priority Queue Is Empty
            while (PQ.Count != 0)
            {
                // Dequeuing The First (Minimum) Element
                // From The PQ
                PQNode Minimum = PQ.Dequeue();

                // Checking If The Minimum Element Is The Root Node Or Not
                // (Only The Root Node Has A Null Parent In The First Iteration)
                if (Minimum.Parent != null)
                {
                    Edge edge = new Edge(Minimum.currentVertex, (RGBPixel)Minimum.Parent, Minimum.currentWeight);
                    MST.Add(edge);     // Add the minimum weight to the MST.
                    MSTSum += edge.weight; // Add The Edge's Weight To The MST Sum
                }

                // We Have To Modify The Values Of The PQ
                // Each Time
                foreach (var node in PQ)
                {
                    // Calculating The Weight Difference Between The Current Node
                    // And The Minimum Node
                    weight = WeightDifference(node, Minimum);

                    // If That Weight Difference Is Less Than The Node's Current Weight
                    // Then The Minimum Node Becomes The Parent Of That Node
                    // And We Adjust That Node's Weight To The Weight Difference
                    // By Updating The PQ
                    if (weight < node.currentWeight)
                    {
                        node.Parent = Minimum.currentVertex;
                        PQ.UpdatePriority(node, weight);
                    }
                }
            }

        }

        #endregion

        #region Clustering Methods
        private static void DFS(int cur, ref HashSet<int> visited, ref Dictionary<int, List<int>> neighbours, ref HashSet<int> cluster)
        {
            visited.Add(cur);
            cluster.Add(cur);
            foreach (var neighbour in neighbours[cur])
            {
                if (!visited.Contains(neighbour))
                    DFS(neighbour, ref visited, ref neighbours, ref cluster);
            }
        }
        public static List<HashSet<int>> Cluster(List<edgee> MST, int k)
        {
            List<HashSet<int>> clusters = new List<HashSet<int>>();
            HashSet<int> visited = new HashSet<int>();
            float MaxWeight;
            int MaxIndex;
            //MST.Sort();
            //int m = k - 1;
            //int size = MST.Count;
            ////foreach (var edge in MST) {
            ////    MessageBox.Show(""+edge.weight);
            ////}
            //for (int i = size-1; i >= m ; i--)
            //{

            //    edgee e = new edgee(MST[i].src,MST[i].dest,MST[i].weight,0);

  
            //    MST[i] = e;
            //}

            for (int j = 0; j < k - 1; j++)
            {
                MaxWeight = 0;
                MaxIndex = 0;

                for (int i = 0; i < MST.Count; i++)
                {
                    if (MST[i].weight > MaxWeight)
                    {
                        MaxWeight = MST[i].weight;
                        MaxIndex = i;
                    }
                }

                edgee e= new edgee(MST[MaxIndex].src, MST[MaxIndex].dest,0);
              
                MST[MaxIndex] = e;
            }

            //Dictionary<RGBPixel, List<RGBPixel>> neighbours = new Dictionary<RGBPixel, List<RGBPixel>>();

            foreach (var edge in nMST)
            {
                if (edge.weight != 0)
                {
                    if (neighbours.ContainsKey(edge.src))
                    {
                        neighbours[edge.src].Add(edge.dest);
                    }
                    else
                    {
                        List<int> l = new List<int>();
                        l.Add(edge.dest);
                        neighbours.Add(edge.src, l);
                    }
                    if (neighbours.ContainsKey(edge.dest))
                    {
                        neighbours[edge.dest].Add(edge.src);
                    }
                    else
                    {
                        List<int> l = new List<int>();
                        l.Add(edge.src);
                        neighbours.Add(edge.dest, l);
                    }
                }
                else
                {
                    if (!neighbours.ContainsKey(edge.src))
                    {
                        List<int> l = new List<int>();
                        neighbours.Add(edge.src, l);
                    }
                    if (!neighbours.ContainsKey(edge.dest))
                    {
                        List<int> l = new List<int>();
                        neighbours.Add(edge.dest, l);
                    }
                }
            }

            int q = 0;

            foreach (var vertex in neighbours)
            {
                if (!visited.Contains(vertex.Key))
                {
                    HashSet<int> h = new HashSet<int>();
                    DFS(vertex.Key, ref visited, ref neighbours, ref h);
                    clusters.Add(h);
                    q++;
                }
            }

            return clusters;
        }

        #endregion

        #region Extracting Representative Colors

        public static List<RGBPixel> ExtractColors(List<HashSet<RGBPixel>> clusters)
        {
            List<RGBPixel> Palette = new List<RGBPixel>();
            RGBPixel avgColor;
            int redSum, greenSum, blueSum;

            for (int i = 0; i < clusters.Count; i++)
            {
                redSum = 0; greenSum = 0; blueSum = 0;

                foreach (var color in clusters[i])
                {
                    redSum += color.red;
                    greenSum += color.green;
                    blueSum += color.blue;
                }

                avgColor.red = (byte)(redSum / clusters[i].Count);
                avgColor.green = (byte)(greenSum / clusters[i].Count);
                avgColor.blue = (byte)(blueSum / clusters[i].Count);

                Palette.Add(avgColor);
            }

            return Palette;
        }

        #endregion

        #region Quantizing The Image

        public static RGBPixel[,] QuantizedImage(RGBPixel[,] imageMatrix, List<RGBPixel> palette)
        {
            int height = imageMatrix.GetLength(0);
            int width = imageMatrix.GetLength(1);

            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    imageMatrix[i, j] = calculateWeight(palette, imageMatrix[i, j]);
                }
            }


            return imageMatrix;
        }
        public static RGBPixel calculateWeight(List<RGBPixel> palette, RGBPixel myCurrentPixel)
        {
            RGBPixel returnPixel = myCurrentPixel;
            double weight = 10000000;
            int r1, g1, b1;
            int r2, g2, b2;

            for (int i = 0; i < palette.Count; i++)
            {
                r1 = palette[i].red;
                g1 = palette[i].green;
                b1 = palette[i].blue;
                r2 = myCurrentPixel.red;
                g2 = myCurrentPixel.green;
                b2 = myCurrentPixel.blue;
                double eq = Math.Sqrt((r2 - r1) * (r2 - r1) + (g2 - g1) * (g2 - g1) + (b2 - b1) * (b2 - b1));

                if (eq < weight)
                {
                    returnPixel = palette[i];
                    weight = eq;
                }
            }
            return returnPixel;
        }

        #endregion
    }
}
