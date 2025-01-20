using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerraMesh.Utils
{
    internal struct EdgePair : IEquatable<EdgePair>
    {
        public readonly int P1;
        public readonly int P2;

        public EdgePair(int p1, int p2)
        {
            P1 = Mathf.Min(p1, p2);
            P2 = Mathf.Max(p1, p2);
        }

        public bool Equals(EdgePair other)
        {
            return P1 == other.P1 && P2 == other.P2;
        }

        public override bool Equals(object obj)
        {
            return obj is EdgePair other && Equals(other);
        }


        public override int GetHashCode()
        {
            return (P1, P2).GetHashCode(); // Use ValueTuple's GetHashCode for efficiency
        }
    }

    internal class QuadTree
    {
        public Bounds bounds;
        public QuadTree[]? children;
        public bool isLeaf = true;

        public QuadTree(Bounds bounds)
        {
            this.bounds = bounds;
        }

        public bool Contains(Vector3 point)
        {
            return bounds.Contains(point);
        }

        public void Subdivide()
        {
            isLeaf = false;
            children = new QuadTree[4];
            float quarterX = bounds.size.x * 0.25f;
            float quarterZ = bounds.size.z * 0.25f;

            // Create four child quadrants
            for (int i = 0; i < 4; i++)
            {
                Vector3 center = bounds.center + new Vector3(
                    ((i % 2) * 2 - 1) * quarterX,
                    0,
                    ((i / 2) * 2 - 1) * quarterZ
                );

                Bounds childBounds = new Bounds(
                    center,
                    new Vector3(bounds.size.x * 0.5f, bounds.size.y, bounds.size.z * 0.5f)
                );

                children[i] = new QuadTree(childBounds);
            }
        }

        public void Subdivide(Bounds levelBounds, Vector2 stepSize, float minCellStep, float maxCellStep, float falloffSpeed, float maxDistance)
        {
            // Get the point relative to the center
            Vector3 closestPoint = levelBounds.ClosestPoint(bounds.center) - levelBounds.center;
            // Size of the side where the closest point is
            float sideSize = Mathf.Max(Mathf.Abs(closestPoint.x), Mathf.Abs(closestPoint.z)); 
            // Size of the side where the quad is
            Vector3 distanceToCenter =  bounds.center - levelBounds.center;
            float distance = Mathf.Max(Mathf.Abs(distanceToCenter.x), Mathf.Abs(distanceToCenter.z));

            float actualCellStep = minCellStep;
            if (distance > sideSize)
            {
                actualCellStep = Mathf.Lerp(minCellStep, maxCellStep, falloffSpeed*(distance - sideSize) / maxDistance);
            }
            actualCellStep = Mathf.Max(actualCellStep, 1);

            // If the quad is too large for desired step size, subdivide
            if (bounds.size.x > actualCellStep * stepSize.x || bounds.size.z > actualCellStep * stepSize.y)
            {
                Subdivide();
                foreach (var child in children!)
                {
                    child.Subdivide(levelBounds, stepSize, minCellStep, maxCellStep, falloffSpeed, maxDistance);
                }
            }
        }

        public void GetLeafNodes(List<QuadTree> leafNodes)
        {
            if (isLeaf)
            {
                leafNodes.Add(this);
            }
            else
            {
                if (children == null)
                    return;
                    
                foreach (var child in children)
                {
                    child.GetLeafNodes(leafNodes);
                }
            }
        } 
    }

    [Serializable]
    public struct TerraMeshConfig
    {
        // Bounding box for target area
        public Bounds? levelBounds;
        public bool useBounds;
        public bool constrainEdges;
        // Mesh subdivision
        public bool subdivideMesh;
        public float baseEdgeLength;
        //Mesh smoothing
        public bool smoothMesh;
        public int smoothingIterations;
        // UVs
        public bool replaceUvs;
        public bool onlyUVs;
        // Renderer mask
        public uint renderingLayerMask;
        // Terrain conversion
        public int minMeshStep;
        public int maxMeshStep;
        public float falloffSpeed;
        public int targetVertexCount;
        public bool carveHoles;
        public bool refineMesh;
        public bool useMeshCollider;
        public bool copyTrees;
        public bool copyDetail;
        public readonly Shader? terraMeshShader;

        /// <summary>
        /// Creates a new TerraMeshConfig object with the specified parameters.
        /// </summary>
        /// <param name="levelBounds">The bounding box for the target area.</param>
        /// <param name="useBounds">Whether to use the level bounds for mesh density.</param>
        /// <param name="constrainEdges">Whether to constrain the mesh edges to the level bounds.</param>
        /// <param name="subdivideMesh">Whether to subdivide the mesh based on the level bounds.</param>
        /// <param name="baseEdgeLength">The base edge length for mesh subdivision.</param>
        /// <param name="smoothMesh">Whether to smooth the mesh after conversion.</param>
        /// <param name="smoothingIterations">The number of smoothing iterations to apply.</param>
        /// <param name="replaceUvs">Whether to replace the UVs of the mesh.</param>
        /// <param name="onlyUVs">Whether to only generate UVs for the mesh.</param>
        /// <param name="renderingLayerMask">The rendering layer mask for the mesh terrain.</param>
        /// <param name="minMeshStep">The (minimum) mesh step for uniform (adaptive) meshing.</param>
        /// <param name="maxMeshStep">The maximum mesh step for adaptive meshing.</param>
        /// <param name="falloffSpeed">The falloff speed for adaptive mesh density.</param>
        /// <param name="targetVertexCount">The target vertex count for adaptive/uniform meshing.</param>
        /// <param name="carveHoles">Whether to carve holes in the mesh terrain.</param>
        /// <param name="refineMesh">Whether to refine the mesh using Triangle.NET quality mesher</param>
        /// <param name="useMeshCollider">Whether to use a MeshCollider for the mesh terrain.</param>
        /// <param name="copyTrees">Whether to copy trees to the mesh terrain.</param>
        /// <param name="copyDetail">Whether to copy detail objects to the mesh terrain. (basic implementation, bad results)</param>
        /// <param name="terraMeshShader">The shader to use for the mesh terrain.</param>
        /// <remarks>
        /// The default shader is MeshTerrainLit.
        /// </remarks>
        public TerraMeshConfig(Bounds? levelBounds = null,
                                bool useBounds = false,
                                bool constrainEdges = true,
                                bool subdivideMesh = true,
                                float baseEdgeLength = 5f,
                                bool smoothMesh = true,
                                int smoothingIterations = 1,
                                bool replaceUvs = false,
                                bool onlyUVs = false,
                                uint renderingLayerMask = 0,
                                int minMeshStep = 1,
                                int maxMeshStep = 32,
                                float falloffSpeed = 3f,
                                int targetVertexCount = -1,
                                bool carveHoles = true,
                                bool refineMesh = true,
                                bool useMeshCollider = false,
                                bool copyTrees = false,
                                bool copyDetail = false,
                                Shader? terraMeshShader = null)
        {
            this.levelBounds = levelBounds;
            this.useBounds = useBounds;
            this.constrainEdges = constrainEdges;
            this.subdivideMesh = subdivideMesh;
            this.baseEdgeLength = baseEdgeLength;
            this.smoothMesh = smoothMesh;
            this.smoothingIterations = smoothingIterations;
            this.replaceUvs = replaceUvs;
            this.onlyUVs = onlyUVs;
            this.renderingLayerMask = renderingLayerMask;
            this.minMeshStep = minMeshStep;
            this.maxMeshStep = maxMeshStep;
            this.falloffSpeed = falloffSpeed;
            this.targetVertexCount = targetVertexCount;
            this.carveHoles = carveHoles;
            this.refineMesh = refineMesh;
            this.useMeshCollider = useMeshCollider;
            this.copyTrees = copyTrees;
            this.copyDetail = copyDetail;
            this.terraMeshShader = terraMeshShader == null ? TerraMeshPlugin.terraMeshShader : terraMeshShader;

            if (this.terraMeshShader == null)
            {
                Debug.LogError("TerraMesh shader not found. This will cause the mesh terrain to have broken visuals.");
            }
        }
    }

    [Serializable]
    public class MeshifyData
    {
        public List<Vector3>? vertices;
        public List<int>? triangles;
        public List<Vector2>? uvs;
        public List<Vector2>? uvs2;

        public MeshifyData(List<Vector3>? vertices = null,
                            List<int>? triangles = null,
                            List<Vector2>? uvs = null,
                            List<Vector2>? uvs2 = null)
        {
            this.vertices = vertices;
            this.triangles = triangles;
            this.uvs = uvs;
            this.uvs2 = uvs2;
        }
    }


    [Serializable]
    public class MeshifyTerrainData: MeshifyData
    {
        public float[,] heightmapData;
        public bool[,] holesData;
        public int holesResolution;
        public int heightmapResolution;
        public float terrainStepX;
        public float terrainStepZ;
        public float terrainWidth;
        public float terrainLength;
        public float terrainHeight;
        public float uvStepX;
        public float uvStepZ;
        public Vector3 terrainPosition;
        public Bounds terrainBounds;

        public MeshifyTerrainData(Terrain terrain)
        {
            GatherTerrainData(terrain,
                                out heightmapData,
                                out holesData,
                                out holesResolution,
                                out heightmapResolution,
                                out terrainStepX,
                                out terrainStepZ,
                                out terrainWidth,
                                out terrainLength,
                                out terrainHeight,
                                out uvStepX,
                                out uvStepZ,
                                out terrainPosition,
                                out terrainBounds);

            vertices = new List<Vector3>();
            triangles = new List<int>();
            uvs = new List<Vector2>();
        }

        private void GatherTerrainData(Terrain terrain,
                                            out float[,] heightmapData,
                                            out bool[,] holesData,
                                            out int holesResolution,
                                            out int heightmapResolution,
                                            out float terrainStepX,
                                            out float terrainStepZ,
                                            out float terrainWidth,
                                            out float terrainLength,
                                            out float terrainHeight,
                                            out float uvStepX,
                                            out float uvStepZ,
                                            out Vector3 terrainPosition,
                                            out Bounds terrainBounds)
        {
            heightmapData = terrain.terrainData.GetHeights(0, 0, terrain.terrainData.heightmapResolution, terrain.terrainData.heightmapResolution);
            holesData = terrain.terrainData.GetHoles(0, 0, terrain.terrainData.holesResolution, terrain.terrainData.holesResolution);
            holesResolution = terrain.terrainData.holesResolution;
            heightmapResolution = terrain.terrainData.heightmapResolution;

            terrainStepX = terrain.terrainData.size.x / (heightmapResolution - 1);
            terrainStepZ = terrain.terrainData.size.z / (heightmapResolution - 1);

            terrainWidth = terrain.terrainData.size.x;
            terrainLength = terrain.terrainData.size.z;
            terrainHeight = terrain.terrainData.size.y;

            uvStepX = 1.0f / (heightmapResolution - 1);
            uvStepZ = 1.0f / (heightmapResolution - 1);

            terrainPosition = terrain.transform.position;
            terrainBounds = terrain.terrainData.bounds;
        }
    }
}