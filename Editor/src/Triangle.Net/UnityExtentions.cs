using System.Collections;
using System.Collections.Generic;
using System.Linq;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using UnityEngine;

namespace TriangleNet.Unity
{
    public static class UnityExtentions
    {
        public static void Add(this Polygon polygon, List<Vector2> contour, bool isHole = false)
        {
            polygon.Add(new Contour(contour.ToTriangleNetVertices()), isHole);
        }

        public static void Add(this Polygon polygon, Vector2 vertex)
        {
            polygon.Add(new Vertex(vertex.x, vertex.y));
        }

        public static Mesh GenerateUnityMesh(this TriangleNetMesh triangleNetMesh, QualityOptions options = null)
        {
            if (options != null)
            {
                triangleNetMesh.Refine(options);
            }
         
            Mesh mesh = new Mesh();
            var triangleNetVerts = triangleNetMesh.Vertices.ToList();
  
            var triangles = triangleNetMesh.Triangles;
       
            Vector3[] verts = new Vector3[triangleNetVerts.Count];
            int[] trisIndex = new int[triangles.Count * 3];

            for (int i = 0; i < verts.Length; i++)
            {
                verts[i] = (Vector3) triangleNetVerts[i];
            }
            
            int k = 0;
         
            foreach (var triangle in triangles)
            {
                for (int i = 2; i >= 0; i--)
                {
                    trisIndex[k] = triangleNetVerts.IndexOf(triangle.GetVertex(i));
                    k++;
                }
            }

            mesh.vertices = verts;
            mesh.triangles = trisIndex;

            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            return mesh;
        }
        
        public static List<Vertex> ToTriangleNetVertices(this List<Vector2> points)

        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (var vec in points)
            {
                vertices.Add(new Vertex(vec.x, vec.y));
            }

            return vertices;
        }

        public static List<Vertex> ToTriangleNetVertices(this List<Vector3> points, int upAxis = 2)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (var vec in points)
            {
                Vertex vertex;
                if (upAxis == 0)
                {
                    vertex = new Vertex(vec.y, vec.z); // X up
#if USE_Z || UNITY_EDITOR
                    vertex.Z = vec.x;
#endif
                }
                else if (upAxis == 1)
                {
                    vertex = new Vertex(vec.x, vec.z); // Y up
#if USE_Z || UNITY_EDITOR
                    vertex.Z = vec.y;
#endif
                }
                else if (upAxis == 2)
                {
                    vertex = new Vertex(vec.x, vec.y); // Z up
#if USE_Z || UNITY_EDITOR
                    vertex.Z = vec.z;
#endif
                }    
                else
                    throw new System.ArgumentException("Invalid upAxis value. Must be 0, 1, or 2.");
                vertices.Add(vertex);
            }

            return vertices;
        }

        public static Vertex ToTriangleNetVertex(this Vector3 point, int upAxis = 2)
        {
            Vertex vertex;
            
            if (upAxis == 0)
            {
                vertex = new Vertex(point.y, point.z); // X up
#if USE_Z || UNITY_EDITOR
                vertex.Z = point.x;
#endif
            }
            else if (upAxis == 1)
            {
                vertex = new Vertex(point.x, point.z); // Y up
#if USE_Z || UNITY_EDITOR
                vertex.Z = point.y;
#endif
            }
            else if (upAxis == 2)
            {
                vertex = new Vertex(point.x, point.y); // Z up
#if USE_Z || UNITY_EDITOR
                vertex.Z = point.z;
#endif
            }    
            else
                throw new System.ArgumentException("Invalid upAxis value. Must be 0, 1, or 2.");

            return vertex;
        }

#if USE_UV || UNITY_EDITOR || UNITY_EDITOR || UNITY_EDITOR
        public static Vertex ToTriangleNetVertex(this Vector3 point, Vector2 uv, int upAxis = 2)
        {
            Vertex vertex;
            
            vertex = point.ToTriangleNetVertex(upAxis);
            vertex.UV = uv;
            return vertex;
        }
#endif

#if USE_Z || UNITY_EDITOR
        public static Vector3 ToVector3(this Vertex vertex, int upAxis = 2)
        {
            Vector3 point;
            if (upAxis == 0)
            {
                point = new Vector3(vertex.Z, vertex.X, vertex.Y); // X up
            }
            else if (upAxis == 1)
            {
                point = new Vector3(vertex.X, vertex.Z, vertex.Y); // Y up
            }
            else if (upAxis == 2)
            {
                point = new Vector3(vertex.X, vertex.Y, vertex.Z); // Z up
            }
            else
                throw new System.ArgumentException("Invalid upAxis value. Must be 0 (X), 1 (Y), or 2 (Z).");
            return point;
        }
#endif

        public static void DrawGizmos(this TriangleNetMesh triangleNetMesh)
        {
            foreach (var triangle in triangleNetMesh.triangles)
            {
                var verts = triangle.vertices;
                Gizmos.DrawLine((Vector3) verts[0], (Vector3) verts[1]);
                Gizmos.DrawLine((Vector3) verts[1], (Vector3) verts[2]);
                Gizmos.DrawLine((Vector3) verts[2], (Vector3) verts[0]);	
            }
        }
    }
}