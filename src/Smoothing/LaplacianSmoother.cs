using TriangleNet.Meshing;
using TriangleNet.Tools;
using TriangleNet.Geometry;
using UnityEngine;
using System.Collections.Generic;

namespace TriangleNet.Smoothing
{
    public class LaplacianSmoother : ISmoother
    {
        float lambda = 0.5f;

        public LaplacianSmoother(float lambda)
        {
            this.lambda = lambda;
        }

        public void Smooth(IMesh mesh)
        {
            Smooth(mesh, 1);
        }

        public void Smooth(IMesh mesh, int limit)
        {
            Smooth(mesh, limit, false);
        }

        // Note: This method is not part of the original Triangle.NET code
        // Not entirely correct since it updates the vertices in place instead of working on a copy
        public void Smooth(IMesh mesh, int limit, bool smoothZ = false)
        {
            var tmesh = mesh as TriangleNetMesh;
            var adj = new AdjacencyMatrix(tmesh);
            var vertices = tmesh.vertices;
            for (int i = 0; i < limit; i++)
            {
                foreach (var vertex in vertices.Values)
                {
                    if (vertex.label == 0) // Only smooth inner vertices
                    {
                        int start = adj.ColumnPointers[vertex.id];
                        int end = adj.ColumnPointers[vertex.id + 1];

                        float _x = 0;
                        float _y = 0;
#if USE_Z
                        float _z = 0;
#endif
#if USE_UV
                        Vector2 _uv = Vector2.zero;
#endif
                        int count = 0;

                        for (int k = start; k < end; k++)
                        {
                            int neighborId = adj.RowIndices[k];
                            if (vertices.TryGetValue(neighborId, out var neighbor))
                            {
                                _x += neighbor.x;
                                _y += neighbor.y;
#if USE_UV
                                _uv += neighbor.UV;
#endif
#if USE_Z
                                if (smoothZ)
                                {
                                    _z += neighbor.z;
                                }
#endif
                                count++;
                            }
                            else
                            {
                                Debug.LogWarning("Neighbor not found");
                            }
                        }

                        if (count > 0)
                        {
                            vertex.x = (1-lambda)*vertex.x +  lambda* _x / count;
                            vertex.y = (1-lambda)*vertex.y +  lambda* _y / count;
#if USE_Z
                            if (smoothZ)
                            {
                                vertex.z = (1-lambda)*vertex.z +  lambda* _z / count;
                            }
#endif
#if USE_UV
                            vertex.UV = (1-lambda)*vertex.UV +  lambda* _uv / count;
#endif
                        }
                        // If no neighbors, coordinates remain unchanged
                    }
                }
            }
        }
    }
}