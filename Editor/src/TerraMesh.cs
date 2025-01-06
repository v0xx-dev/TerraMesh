using System;
using System.Linq;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Rendering;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Unity;
using TriangleNet.Smoothing;

#if UNITY_EDITOR
using UnityEditor;
using System.IO;
#endif

namespace TerraMesh
{
#if UNITY_EDITOR
    public class TerraMesh : EditorWindow
    {
        [Tooltip("The Terrain object to convert to a mesh.")]
        private Terrain terrain;
        [Tooltip("The BoxCollider object representing the level bounds.")]
        private BoxCollider levelBounds;
        [Tooltip("The minimum mesh step size. Smaller values result in higher resolution meshes.")]
        private int minMeshStep = 1;
        [Tooltip("The maximum mesh step size.")]
        private int maxMeshStep = 16;
        [Tooltip("How quickly the mesh step size transitions from min to max step size with distance from the level bounds.")]
        private float falloffSpeed = 2f;
        [Tooltip("Refine the mesh by subdividing thin triangles.")]
        private bool refineMesh = true;
        [Tooltip("Use a mesh collider for the mesh terrain.")]
        private bool useMeshCollider = true;
        [Tooltip("Carve holes in the mesh terrain.")]
        private bool carveHoles = true;
        [Tooltip("Copy trees from the terrain to the mesh terrain.")]
        private bool copyTrees = false;
        private bool copyDetail = false;
        
        private int targetVertexCount = -1;

        

        private string saveFolderPath = "Assets/TerrainMeshes"; // Default save folder
        private GameObject meshTerrain;

        [MenuItem("Tools/TerraMesh")]
        public static void ShowWindow()
        {
            GetWindow<TerraMesh>("TerraMesh");
        }

        private void OnGUI()
        {
            EditorGUILayout.LabelField("Basic Settings", EditorStyles.boldLabel);
            
            EditorGUILayout.BeginHorizontal();
            saveFolderPath = EditorGUILayout.TextField("Save Folder", saveFolderPath);

            if (GUILayout.Button("Browse", GUILayout.Width(60)))
            {
                string selectedPath = EditorUtility.OpenFolderPanel("Select Save Folder", saveFolderPath, "");
                if (!string.IsNullOrEmpty(selectedPath))
                {
                    // Make path relative to the project's Assets folder
                    string projectPath = Application.dataPath;
                    if (selectedPath.StartsWith(projectPath))
                    {
                        saveFolderPath = "Assets" + selectedPath.Substring(projectPath.Length);
                    }
                    else
                    {
                        Debug.LogError("Selected folder must be within the project's Assets folder!");
                    }
                }
            }
            EditorGUILayout.EndHorizontal();

            terrain = (Terrain)EditorGUILayout.ObjectField("Terrain", terrain, typeof(Terrain), true);
            levelBounds = (BoxCollider)EditorGUILayout.ObjectField("Level Bounds", levelBounds, typeof(BoxCollider), true);
            //targetVertexCount = EditorGUILayout.IntField("Target Vertex Count", targetVertexCount);
            EditorGUILayout.LabelField("Sampling Settings", EditorStyles.boldLabel);
            minMeshStep = EditorGUILayout.IntField("Minimum Mesh Step", minMeshStep);
            if (levelBounds != null)
            {
                maxMeshStep = EditorGUILayout.IntField("Maximum Mesh Step", maxMeshStep);
                falloffSpeed = EditorGUILayout.FloatField("Falloff Speed", falloffSpeed);
                refineMesh = EditorGUILayout.Toggle("Refine Mesh", refineMesh);
            }
            EditorGUILayout.LabelField("Copy Settings", EditorStyles.boldLabel);
            useMeshCollider = EditorGUILayout.Toggle("Mesh Collider", useMeshCollider);
            carveHoles = EditorGUILayout.Toggle("Carve Holes", carveHoles);
            copyTrees = EditorGUILayout.Toggle("Copy Trees", copyTrees);
            //copyDetail = EditorGUILayout.Toggle("Copy Detail", copyDetail);

            if (GUILayout.Button("Turn to Mesh"))
            {
                
                if (meshTerrain != null)
                {
                    DestroyImmediate(meshTerrain);
                }

                if (ValidateSettings())
                {
                    Shader? terraMeshShader = Shader.Find("Shader Graphs/MeshTerrainLit");

                    if (terraMeshShader == null)
                    {
                        Debug.LogError("Couldn't find the mesh terrain shader!");
                        return;
                    }

                    string currentTime = DateTime.Now.ToString("ddMHHmmss"); 
                    string terrainFolderPath = Path.Combine(saveFolderPath, terrain.name + "_" + currentTime); // Subfolder for the terrain

                    // Create subfolder if it doesn't exist
                    Directory.CreateDirectory(terrainFolderPath);

                    TerraMeshConfig config = new TerraMeshConfig(
                        levelBounds : levelBounds?.bounds,
                        useMeshCollider : useMeshCollider,
                        targetVertexCount : targetVertexCount,
                        minMeshStep : minMeshStep,
                        maxMeshStep : maxMeshStep,
                        falloffSpeed : falloffSpeed,
                        refineMesh : refineMesh,
                        carveHoles : carveHoles,
                        copyTrees : copyTrees,
                        copyDetail : copyDetail,
                        terraMeshShader : terraMeshShader
                    );

                    meshTerrain = terrain.Meshify(config);

                    string materialPath = Path.Combine(terrainFolderPath, terrain.name + "_Lit.mat");
                    Material terrainMaterial = new Material(terraMeshShader); 
                    AssetDatabase.CreateAsset(terrainMaterial, materialPath);
                    meshTerrain.GetComponent<MeshRenderer>().sharedMaterial = terrainMaterial;
                    
                    terrainMaterial.SetupMaterialFromTerrain(terrain, terrainFolderPath);

                    string meshPath = Path.Combine(terrainFolderPath, terrain.name + "_Mesh.asset");
                    AssetDatabase.CreateAsset(meshTerrain.GetComponent<MeshFilter>().sharedMesh, meshPath);

                    string prefabPath = Path.Combine(terrainFolderPath, terrain.name + ".prefab");
                    PrefabUtility.SaveAsPrefabAsset(meshTerrain, prefabPath);
                }

            }
        }

        private bool ValidateSettings()
        {
            if (terrain == null)
            {
                Debug.LogError("Terrain is not set!");
                return false;
            }
            if (levelBounds == null)
            {
                Debug.LogWarning("Level bounds are not set. Uniform meshing will be applied.");
            }

            if (minMeshStep < 1 || maxMeshStep < 1)
            {
                Debug.LogWarning("Mesh step should be at positive and at least one! Setting to 1 to avoid errors...");
                minMeshStep = Mathf.Max(minMeshStep, 1);
                maxMeshStep = Mathf.Max(maxMeshStep, 1);
            }

            return true;
        }
    }
#endif

    public static class TerraMeshExtensions
    {
        /// <summary>
        /// Postprocesses a mesh terrain by refining, smoothing, and unwrapping its UVs.
        /// </summary>
        /// <param name="meshTerrainObject">The GameObject containing the mesh terrain to postprocess.</param>
        /// <param name="config">The TerraMeshConfig object containing the postprocessing parameters.</param>
        /// <remarks>
        /// This method refines the mesh terrain by subdividing triangles longer than the specified edge length.
        /// It then smooths the mesh terrain using Laplacian smoothing with the specified number of iterations.
        /// Finally, it unwraps the UVs of the mesh terrain using the specified height axis.
        /// </remarks>
        /// <exception cref="System.ArgumentNullException">Thrown if `meshTerrainObject` or `config` is null.</exception>
        /// <exception cref="System.ArgumentException">Thrown if `config` contains invalid parameters.</exception>
        public static void PostprocessMeshTerrain(this GameObject meshTerrainObject, TerraMeshConfig config)
        {
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();

            MeshFilter? meshFilter = meshTerrainObject!.GetComponent<MeshFilter>();
            MeshCollider? meshCollider = meshTerrainObject.GetComponent<MeshCollider>();
            MeshRenderer? meshRenderer = meshTerrainObject.GetComponent<MeshRenderer>();
            
            Bounds levelBounds = config.levelBounds.HasValue ? config.levelBounds.Value : meshRenderer?.bounds ?? default;

            if (meshFilter == null || meshCollider == null)
            {
                Debug.LogError("MeshFilter and MeshCollider components must be present on the object to modify.");
            }

            if (levelBounds == default)
            {
                Debug.LogError("Level bounds must be set to postprocess the mesh terrain. Either set the level bounds in the TerraMeshConfig or ensure the MeshRenderer component is present on the object.");
                return;
            }

            Mesh originalMesh = meshFilter!.sharedMesh;
            (Mesh newMesh, int submeshIndex) = meshTerrainObject.ExtractLargestSubmesh();
            newMesh.name = originalMesh.name + "_Postprocessed";

            Debug.LogDebug("Extracted submesh with " + newMesh.vertexCount + " vertices and " + newMesh.triangles.Length / 3 + " triangles" + " from submesh " + submeshIndex + " of " + meshTerrainObject.name);

            // Determine height axis
            int heightAxis = 0;

            Vector3 objectSpaceUpDirection = meshTerrainObject.transform.InverseTransformDirection(Vector3.up);
            objectSpaceUpDirection.x *= Mathf.Sign(meshTerrainObject.transform.lossyScale.x);
            objectSpaceUpDirection.y *= Mathf.Sign(meshTerrainObject.transform.lossyScale.y);
            objectSpaceUpDirection.z *= Mathf.Sign(meshTerrainObject.transform.lossyScale.z);

            if (Mathf.Abs(Vector3.Dot(objectSpaceUpDirection, Vector3.right)) > 0.5f)
            {
                heightAxis = 0;
            }
            else if (Mathf.Abs(Vector3.Dot(objectSpaceUpDirection, Vector3.up)) > 0.5f)
            {
                heightAxis = 1;
            }
            else if (Mathf.Abs(Vector3.Dot(objectSpaceUpDirection, Vector3.forward)) > 0.5f)
            {
                heightAxis = 2;
            }

            Debug.LogDebug("Object's up axis: " + heightAxis + ". Object's up direction: " + objectSpaceUpDirection);

            if (config.onlyUVs)
            {
                newMesh.uv2 = UnwrapUVs(newMesh.vertices, heightAxis, true);
                meshFilter = meshTerrainObject.GetComponent<MeshFilter>();
                meshFilter.sharedMesh = newMesh;

                meshCollider = meshTerrainObject.GetComponent<MeshCollider>();
                meshCollider.sharedMesh = newMesh;

                return;
            }

            Bounds objectSpaceBounds;
            if (config.useBounds)
            {
                // Transform bounds to object space
                objectSpaceBounds = new Bounds(
                    meshTerrainObject.transform.InverseTransformPoint(levelBounds.center),
                    meshTerrainObject.transform.InverseTransformVector(levelBounds.size)
                );
                //Ensure bounds are not negative (take absolute value)
                objectSpaceBounds.size = new Vector3(Mathf.Abs(objectSpaceBounds.size.x), Mathf.Abs(objectSpaceBounds.size.y), Mathf.Abs(objectSpaceBounds.size.z));
                //Check if object space bounds intersect with mesh bounds
                if (!newMesh.bounds.Intersects(objectSpaceBounds))
                {
                    Debug.LogDebug("Object's bounds do not intersect with level bounds. No vertices will be processed.");
                    return;
                }
            }
            else
            {
                objectSpaceBounds = newMesh.bounds;
            }

            Debug.LogDebug("Bounds for triangulation: " + objectSpaceBounds.center + " " + objectSpaceBounds.size);

            // Get the original mesh data
            Vector3[] vertices = newMesh.vertices;
            Vector2[] uvs = newMesh.uv;
            Vector2[] uvs2 = new Vector2[0];
            int[] triangles = newMesh.triangles;

            //Calculate max scale only taking the absolute horizontal scales into account, horizontal is everything not on the height axis
            float maxScale = Mathf.Max(Mathf.Abs(meshTerrainObject.transform.lossyScale[(heightAxis + 1) % 3]), Mathf.Abs(meshTerrainObject.transform.lossyScale[(heightAxis + 2) % 3]));
            float maxEdgeLength = config.baseEdgeLength / maxScale;
            Debug.LogDebug("Max scale: " + maxScale + " Max edge length: " + maxEdgeLength + "m");

            BeautifyMeshData(config, ref vertices, ref triangles, ref uvs, ref uvs2, maxEdgeLength, heightAxis, objectSpaceBounds, sw);

            newMesh.vertices = vertices;
            newMesh.triangles = triangles;
            newMesh.uv = uvs;
            newMesh.uv2 = uvs2;

            newMesh.Optimize();
            newMesh.RecalculateNormals();
            newMesh.RecalculateTangents();
            newMesh.RecalculateBounds();

            meshFilter.sharedMesh = newMesh;
            meshCollider!.sharedMesh = newMesh;
        }

        private static void BeautifyMeshData(
                                            TerraMeshConfig config,
                                            ref Vector3[] vertices,
                                            ref int[] triangles,
                                            ref Vector2[] uvs,
                                            ref Vector2[] uvs2,
                                            float maxEdgeLength,
                                            int heightAxis,
                                            Bounds objectSpaceBounds,
                                            System.Diagnostics.Stopwatch sw)
        {
            // Create lists to store vertices and triangles within bounds
            HashSet<int> vertexIndicesInBounds = new HashSet<int>();
            Queue<(int, int, int)> trianglesToRefine = new Queue<(int, int, int)>();
            // Track unique vertices to avoid duplicates
            Dictionary<Vector3, int> uniqueVertices = new Dictionary<Vector3, int>();
            //Triangles that are outside bounds
            List<int> trianglesToKeep = new List<int>();
            //Edges within bounds
            Dictionary<EdgePair, int> edgeCounts = new Dictionary<EdgePair, int>();
            // Dictionaries to map Triangle.NET vertices to original mesh indices and vice versa
            Dictionary<int, int> poly2VertMap = new Dictionary<int, int>();
            Dictionary<int, int> vert2PolyMap = new Dictionary<int, int>();

            // Find vertices within bounds and triangles to keep
            // Keep triangles that are outside bounds
            Debug.LogDebug("Filtering vertices and triangles... Count: " + triangles.Length / 3 + " Vertices: " + vertices.Length);
            sw.Start();
            for (int i = 0; i < triangles.Length; i += 3)
            {
                int p1 = triangles[i];
                int p2 = triangles[i + 1];
                int p3 = triangles[i + 2];

                //find unique vertices and if found a duplicate vertex, get a reference to the original vertex
                (p1, p2, p3) = CheckForDuplicateVertices(uniqueVertices, vertices, p1, p2, p3);

                bool allVerticesInBounds = objectSpaceBounds.Contains(vertices[p1]) &&
                                            objectSpaceBounds.Contains(vertices[p2]) &&
                                            objectSpaceBounds.Contains(vertices[p3]);

                if (allVerticesInBounds)
                {

                    // Store vertices within bounds
                    vertexIndicesInBounds.Add(p1);
                    vertexIndicesInBounds.Add(p2);
                    vertexIndicesInBounds.Add(p3);

                    // Store triangles within bounds
                    if (config.subdivideMesh)
                    {
                        trianglesToRefine.Enqueue((p1, p2, p3));
                    }

                    // Count edges
                    UpdateEdgeCounts(edgeCounts, p1, p2);
                    UpdateEdgeCounts(edgeCounts, p2, p3);
                    UpdateEdgeCounts(edgeCounts, p3, p1);
                }
                else
                {
                    trianglesToKeep.Add(p1);
                    trianglesToKeep.Add(p2);
                    trianglesToKeep.Add(p3);
                }
            }
            sw.Stop();
            Debug.LogDebug("Vertex and triangle filtering/deduplication time: " + sw.ElapsedMilliseconds + "ms");
            Debug.LogDebug("Vertices in bounds: " + vertexIndicesInBounds.Count + " Triangles in bounds: " + trianglesToRefine.Count);

            // Allow addition of new vertices
            List<Vector3> newVertices = new List<Vector3>(vertices);
            List<Vector2> newUVs = new List<Vector2>(uvs);

            // Refine triangles
            if (config.subdivideMesh)
            {
                Debug.LogDebug($"Refining {trianglesToRefine.Count} triangles...");
                sw.Restart();
                Dictionary<Vector3, int> midpointVertices = new Dictionary<Vector3, int>();
                while (trianglesToRefine.Count > 0)
                {
                    (int p1, int p2, int p3) = trianglesToRefine.Dequeue();

                    Vector3 v1 = newVertices[p1];
                    Vector3 v2 = newVertices[p2];
                    Vector3 v3 = newVertices[p3];

                    // Calculate edge lengths
                    float edge1Length2 = (v1 - v2).sqrMagnitude;
                    float edge2Length2 = (v2 - v3).sqrMagnitude;
                    float edge3Length2 = (v3 - v1).sqrMagnitude;

                    EdgePair edge1 = new EdgePair(p1, p2);
                    EdgePair edge2 = new EdgePair(p2, p3);
                    EdgePair edge3 = new EdgePair(p3, p1);

                    float maxEdgeLength2 = maxEdgeLength * maxEdgeLength;

                    // Check if all of the edges are longer than the maximum length and if none of them are on the boundary
                    bool refine = edge1Length2 > maxEdgeLength2 && edge2Length2 > maxEdgeLength2 && edge3Length2 > maxEdgeLength2;
                    bool onBoundary = edgeCounts.ContainsKey(edge1) && edgeCounts[edge1] == 1 ||
                                    edgeCounts.ContainsKey(edge2) && edgeCounts[edge2] == 1 ||
                                    edgeCounts.ContainsKey(edge3) && edgeCounts[edge3] == 1;

                    if (refine && !onBoundary)
                    {
                        Vector3 midpoint1 = (v1 + v2) / 2;
                        Vector3 midpoint2 = (v2 + v3) / 2;
                        Vector3 midpoint3 = (v3 + v1) / 2;

                        Vector2 uv1 = newUVs[p1];
                        Vector2 uv2 = newUVs[p2];
                        Vector2 uv3 = newUVs[p3];

                        // Add new vertices while keeping track of duplicates
                        int newVertInd1 = newVertices.Count;
                        if (!midpointVertices.ContainsKey(midpoint1))
                        {
                            midpointVertices[midpoint1] = newVertInd1;
                            vertexIndicesInBounds.Add(newVertInd1);
                            newVertices.Add(midpoint1);
                            newUVs.Add((uv1 + uv2) / 2);
                        }
                        else
                        {
                            newVertInd1 = midpointVertices[midpoint1];
                        }

                        int newVertInd2 = newVertices.Count;
                        if (!midpointVertices.ContainsKey(midpoint2))
                        {
                            midpointVertices[midpoint2] = newVertInd2;
                            vertexIndicesInBounds.Add(newVertInd2);
                            newVertices.Add(midpoint2);
                            newUVs.Add((uv2 + uv3) / 2);
                        }
                        else
                        {
                            newVertInd2 = midpointVertices[midpoint2];
                        }

                        int newVertInd3 = newVertices.Count;
                        if (!midpointVertices.ContainsKey(midpoint3))
                        {
                            midpointVertices[midpoint3] = newVertInd3;
                            vertexIndicesInBounds.Add(newVertInd3);
                            newVertices.Add(midpoint3);
                            newUVs.Add((uv3 + uv1) / 2);
                        }
                        else
                        {
                            newVertInd3 = midpointVertices[midpoint3];
                        }

                        // Add new triangles
                        trianglesToRefine.Enqueue((p1, newVertInd1, newVertInd3));
                        trianglesToRefine.Enqueue((newVertInd1, p2, newVertInd2));
                        trianglesToRefine.Enqueue((newVertInd2, p3, newVertInd3));
                        trianglesToRefine.Enqueue((newVertInd1, newVertInd2, newVertInd3));
                    }
                }
                sw.Stop();
                Debug.LogDebug("Refinement time: " + sw.ElapsedMilliseconds + "ms");
            }

            // Create Triangle.NET polygon input
            var polygon = new Polygon();
            // Add vertices within bounds to the polygon
            foreach (int i in vertexIndicesInBounds)
            {
                Vertex triNetVertex = newVertices[i].ToTriangleNetVertex(newUVs[i], heightAxis);
                polygon.Add(triNetVertex);
                poly2VertMap[polygon.Points.Count - 1] = i;
                vert2PolyMap[i] = polygon.Points.Count - 1;
            }

            // Add segments to enforce edges on the outer boundary
            if (config.constrainEdges)
            {
                foreach (var edge in edgeCounts)
                {
                    if (edge.Value == 1)
                    {
                        polygon.Add(new Segment(polygon.Points[vert2PolyMap[edge.Key.P1]], polygon.Points[vert2PolyMap[edge.Key.P2]]));
                    }
                }
            }

            // Configure triangulation options
            int maxAdditionalVertices = Mathf.Min(vertexIndicesInBounds.Count / 4, 30000);
            ConstraintOptions options = new ConstraintOptions() { ConformingDelaunay = true, SegmentSplitting = 2 };
            QualityOptions quality = new QualityOptions() { MinimumAngle = 30.0f, SteinerPoints = config.subdivideMesh ? maxAdditionalVertices : 0 };

            // Perform triangulation
            Debug.LogDebug($"Triangulating {polygon.Points.Count} vertices...");
            sw.Restart();
            var mesh2d = polygon.Triangulate(options, quality);
            sw.Stop();
            Debug.LogDebug("Triangulation time: " + sw.ElapsedMilliseconds + "ms");
            if (config.smoothMesh)
            {
                sw.Restart();
                var smoother = new LaplacianSmoother(1f);
                smoother.Smooth(mesh2d, config.smoothingIterations);
                sw.Stop();
                Debug.LogDebug("Smoothing time: " + sw.ElapsedMilliseconds + "ms");
            }

            // Create new triangles list combining kept triangles and new triangulation
            List<int> newTriangles = new List<int>(trianglesToKeep);

            Debug.LogDebug("Modified triangles: " + mesh2d.Triangles.Count);

            // Convert Triangle.NET triangles back to Unity mesh triangles
            foreach (ITriangle triangle in mesh2d.Triangles)
            {

                int v0 = triangle.GetVertexID(0);
                int v1 = triangle.GetVertexID(1);
                int v2 = triangle.GetVertexID(2);

                // Use the Dictionary to get the original mesh indices
                // If the vertex is new add it to the list

                if (!poly2VertMap.ContainsKey(v0))
                {
                    poly2VertMap[v0] = newVertices.Count;
                    newVertices.Add(triangle.GetVertex(0).ToVector3(heightAxis));
                    newUVs.Add(triangle.GetVertex(0).UV);
                }
                if (!poly2VertMap.ContainsKey(v1))
                {
                    poly2VertMap[v1] = newVertices.Count;
                    newVertices.Add(triangle.GetVertex(1).ToVector3(heightAxis));
                    newUVs.Add(triangle.GetVertex(1).UV);
                }
                if (!poly2VertMap.ContainsKey(v2))
                {
                    poly2VertMap[v2] = newVertices.Count;
                    newVertices.Add(triangle.GetVertex(2).ToVector3(heightAxis));
                    newUVs.Add(triangle.GetVertex(2).UV);
                }

                int vertexIndex0 = poly2VertMap[v0];
                int vertexIndex1 = poly2VertMap[v1];
                int vertexIndex2 = poly2VertMap[v2];

                //Different winding order based on a height axis
                if (heightAxis == 1)
                {
                    newTriangles.Add(vertexIndex0);
                    newTriangles.Add(vertexIndex2);
                    newTriangles.Add(vertexIndex1);
                }
                else
                {
                    newTriangles.Add(vertexIndex0);
                    newTriangles.Add(vertexIndex1);
                    newTriangles.Add(vertexIndex2);
                }

            }

            Debug.LogDebug("Original UVs: " + uvs.Length + " New UVs: " + newUVs.Count);
            Debug.LogDebug("Original vertices: " + vertices.Length + " New vertices: " + newVertices.Count);

            vertices = newVertices.ToArray();
            triangles = newTriangles.ToArray();
            uvs2 = UnwrapUVs(vertices, heightAxis, true);
            
            if (config.replaceUvs)
            {
                uvs = uvs2.ToArray(); // Replace UVs with unwrapped UVs (only acceptable on moons with no splatmaps)
            }
            else
            {
                uvs = newUVs.ToArray();
            }
        }

        private static Vector2[] UnwrapUVs(Vector3[] vertices, int upAxis, bool normalize)
        {
            Vector2[] uvs = new Vector2[vertices.Length];
            Vector2 min = new Vector2(float.MaxValue, float.MaxValue);
            Vector2 max = new Vector2(float.MinValue, float.MinValue);

            for (int i = 0; i < vertices.Length; i++)
            {
                if (upAxis == 0)
                    uvs[i] = new Vector2(vertices[i].y, vertices[i].z);
                else if (upAxis == 1)
                    uvs[i] = new Vector2(vertices[i].x, vertices[i].z);
                else
                    uvs[i] = new Vector2(vertices[i].x, vertices[i].y);

                min = Vector2.Min(min, uvs[i]);
                max = Vector2.Max(max, uvs[i]);
            }
            if (normalize)
            {
                Vector2 size = max - min;
                for (int i = 0; i < uvs.Length; i++)
                {
                    uvs[i] = new Vector2((uvs[i].x - min.x) / size.x, (uvs[i].y - min.y) / size.y);
                }
            }

            return uvs;

        }

        private static void UpdateEdgeCounts(Dictionary<EdgePair, int> edgeCounts, int p1, int p2)
        {
            EdgePair edge = new EdgePair(p1, p2);
            if (edgeCounts.ContainsKey(edge))
                edgeCounts[edge]++;
            else
                edgeCounts[edge] = 1;
        }

        private static (int, int, int) CheckForDuplicateVertices(Dictionary<Vector3, int> uniqueVertices, Vector3[] vertices, int p1, int p2, int p3)
        {
            if (!uniqueVertices.ContainsKey(vertices[p1]))
            {
                uniqueVertices[vertices[p1]] = p1;
            }
            else
            {
                p1 = uniqueVertices[vertices[p1]];
            }
            if (!uniqueVertices.ContainsKey(vertices[p2]))
            {
                uniqueVertices[vertices[p2]] = p2;
            }
            else
            {
                p2 = uniqueVertices[vertices[p2]];
            }
            if (!uniqueVertices.ContainsKey(vertices[p3]))
            {
                uniqueVertices[vertices[p3]] = p3;
            }
            else
            {
                p3 = uniqueVertices[vertices[p3]];
            }

            return (p1, p2, p3);
        }

        public static (Mesh submesh, int submeshIndex) ExtractLargestSubmesh(this GameObject meshObject)
        {
            Mesh mesh = meshObject.GetComponent<MeshFilter>().sharedMesh;
            Transform transform = meshObject.transform;

            Mesh meshCopy = mesh.MakeReadableCopy();

            var submeshCount = meshCopy.subMeshCount;
            if (submeshCount <= 1) return (meshCopy, 0);


            var largestSubmeshIndex = Enumerable.Range(0, submeshCount)
                .OrderByDescending(i => meshCopy.GetSubMesh(i).vertexCount)
                .First();

            var triangles = meshCopy.GetTriangles(largestSubmeshIndex);
            var usedVertices = new HashSet<int>(triangles);
            var oldToNewVertexMap = new Dictionary<int, int>();

            Vector3[] vertices = meshCopy.vertices;
            Vector2[] uv = meshCopy.uv;
            Vector3[] normals = meshCopy.normals;
            Vector4[] tangents = meshCopy.tangents;

            var newVertices = new List<Vector3>();
            var newNormals = meshCopy.normals.Length > 0 ? new List<Vector3>() : null;
            var newUVs = meshCopy.uv.Length > 0 ? new List<Vector2>() : null;
            var newTangents = meshCopy.tangents.Length > 0 ? new List<Vector4>() : null;

            foreach (var oldIndex in usedVertices)
            {
                oldToNewVertexMap[oldIndex] = newVertices.Count;
                newVertices.Add(transform.InverseTransformPoint(vertices[oldIndex]));
                if (newNormals != null) newNormals.Add(normals[oldIndex]);
                if (newUVs != null) newUVs.Add(uv[oldIndex]);
                if (newTangents != null) newTangents.Add(tangents[oldIndex]);
            }

            var newTriangles = triangles.Select(oldIndex => oldToNewVertexMap[oldIndex]).ToArray();

            var submesh = new Mesh
            {
                name = meshCopy.name + "_Submesh" + largestSubmeshIndex,
                vertices = newVertices.ToArray(),
                triangles = newTriangles,
                indexFormat = meshCopy.indexFormat
            };

            if (newNormals != null) submesh.normals = newNormals.ToArray();
            if (newUVs != null) submesh.uv = newUVs.ToArray();
            if (newTangents != null) submesh.tangents = newTangents.ToArray();

            submesh.RecalculateBounds();
#if UNITY_EDITOR
            GameObject.DestroyImmediate(meshCopy);
#else
            GameObject.Destroy(meshCopy);
#endif

            return (submesh, largestSubmeshIndex);
        }

        //Credit to Matty for this method
        public static Mesh MakeReadableCopy(this Mesh nonReadableMesh)
        {
            if (nonReadableMesh.isReadable)
                return nonReadableMesh;

            var meshCopy = new Mesh();
            meshCopy.indexFormat = nonReadableMesh.indexFormat;

            // Handle vertices
            nonReadableMesh.vertexBufferTarget = GraphicsBuffer.Target.Vertex;
            if (nonReadableMesh.vertexBufferCount > 0)
            {
                var verticesBuffer = nonReadableMesh.GetVertexBuffer(0);
                var totalSize = verticesBuffer.stride * verticesBuffer.count;
                var data = new byte[totalSize];
                verticesBuffer.GetData(data);
                meshCopy.SetVertexBufferParams(nonReadableMesh.vertexCount, nonReadableMesh.GetVertexAttributes());
                meshCopy.SetVertexBufferData(data, 0, 0, totalSize);
                verticesBuffer.Release();
            }

            // Handle triangles
            nonReadableMesh.indexBufferTarget = GraphicsBuffer.Target.Index;
            meshCopy.subMeshCount = nonReadableMesh.subMeshCount;
            var indexesBuffer = nonReadableMesh.GetIndexBuffer();
            var tot = indexesBuffer.stride * indexesBuffer.count;
            var indexesData = new byte[tot];
            indexesBuffer.GetData(indexesData);
            meshCopy.SetIndexBufferParams(indexesBuffer.count, nonReadableMesh.indexFormat);
            meshCopy.SetIndexBufferData(indexesData, 0, 0, tot);
            indexesBuffer.Release();

            // Restore submesh structure
            uint currentIndexOffset = 0;
            for (var i = 0; i < meshCopy.subMeshCount; i++)
            {
                var subMeshIndexCount = nonReadableMesh.GetIndexCount(i);
                meshCopy.SetSubMesh(i, new SubMeshDescriptor((int)currentIndexOffset, (int)subMeshIndexCount));
                currentIndexOffset += subMeshIndexCount;
            }

            // Recalculate normals and bounds
            meshCopy.RecalculateNormals();
            meshCopy.RecalculateBounds();

            meshCopy.name = $"Readable {nonReadableMesh.name}";
            return meshCopy;
        }

        /// <summary>
        /// Converts a Terrain object to a mesh using the specified configuration.
        /// </summary>
        /// <param name="terrain">The Terrain object to convert to a mesh.</param>
        /// <param name="config">The TerraMeshConfig object containing the mesh conversion parameters.</param>
        /// <returns>The GameObject containing the mesh terrain.</returns>
        /// <remarks>
        /// This method creates a new GameObject with a MeshFilter and MeshRenderer component, and assigns the mesh terrain to the MeshFilter.
        /// The mesh terrain is created by sampling the terrain heightmap and converting it to a mesh.
        /// The mesh terrain is parented to the same GameObject as the Terrain object, and has the same layer, tag, and rendering layer mask.
        /// </remarks>
        /// <exception cref="System.ArgumentNullException">Thrown if `terrain` is null.</exception>
        public static GameObject Meshify(this Terrain terrain, TerraMeshConfig config)
        {
            if (terrain == null)
            {
                throw new ArgumentNullException("terrain", "Terrain object cannot be null.");
            }

            MeshifyTerrainData meshTerrainData = new MeshifyTerrainData(terrain);

            GenerateMeshData(config, meshTerrainData);

            GameObject meshTerrain = terrain.Meshify(config, meshTerrainData);

            return meshTerrain;
        }

        private static GameObject Meshify(this Terrain terrain, TerraMeshConfig config, MeshifyTerrainData meshTerrainData)
        {
            Vector3[] vertices = meshTerrainData.vertices.ToArray();
            int[] triangles = meshTerrainData.triangles.ToArray(); 
            Vector2[] uvs = meshTerrainData.uvs.ToArray();

            // Create mesh from sampled vertices and triangles
            Mesh mesh = new Mesh();
            mesh.name = "MeshTerrain_" + terrain.name;
            mesh.indexFormat = vertices.Length > 65535 ? UnityEngine.Rendering.IndexFormat.UInt32 : UnityEngine.Rendering.IndexFormat.UInt16;
            mesh.vertices = vertices;
            mesh.uv = uvs;
            mesh.triangles = triangles;
            mesh.Optimize(); // Clean unused vertices
            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            mesh.RecalculateTangents();

            // Create GameObject with MeshFilter and MeshRenderer components
            GameObject meshTerrain = new GameObject("MeshTerrain_" + terrain.name);
            MeshFilter meshFilter = meshTerrain.AddComponent<MeshFilter>();
            MeshRenderer meshRenderer = meshTerrain.AddComponent<MeshRenderer>();
            // Set same position, parent, rendering layer, rendering layer mask and tag as the terrain (and set the snow overlay custom pass layer)
            meshTerrain.transform.position = terrain.transform.position;
            meshTerrain.transform.SetParent(terrain.transform.parent);
            meshRenderer.gameObject.layer = terrain.gameObject.layer;
            meshRenderer.gameObject.tag = terrain.gameObject.tag;
            terrain.renderingLayerMask |= config.renderingLayerMask;
            meshRenderer.renderingLayerMask = terrain.renderingLayerMask;
            meshTerrain.isStatic = true;

            if (config.useMeshCollider)
            {
                MeshCollider meshCollider = meshTerrain.AddComponent<MeshCollider>();
                meshCollider.sharedMesh = mesh;
                // Disable terrain collider
                terrain.GetComponent<TerrainCollider>().enabled = false;
            }

            meshFilter.mesh = mesh;
            meshRenderer.sharedMaterial = new Material(config.terraMeshShader);
            // Disable rendering of terrain
            terrain.drawHeightmap = false;
            // terrain.drawInstanced = false; // TODO: Check if this is necessary

            //Ideally trees should be copied if we want to use a mesh collider
            if (config.copyTrees)
            {
                //Disable rendering of trees on terrain
                terrain.treeDistance = 0;
                // Create Trees parent object
                Transform treesParent = new GameObject("Trees").transform;
                treesParent.position = meshTerrain.transform.position;
                treesParent.parent = meshTerrain.transform;

                // Copy trees to the mesh terrain
                foreach (TreeInstance tree in terrain.terrainData.treeInstances)
                {
                    // Get tree prototype
                    GameObject treePrototype = terrain.terrainData.treePrototypes[tree.prototypeIndex].prefab;

                    // Instantiate tree
                    GameObject newTree = GameObject.Instantiate(treePrototype, treesParent);

                    // Calculate tree terrainPosition
                    Vector3 treeWorldPos = Vector3.Scale(tree.position, terrain.terrainData.size) + terrain.transform.position;

                    // Set tree transform
                    newTree.transform.position = treeWorldPos;
                    newTree.transform.localScale = new Vector3(tree.widthScale, tree.heightScale, tree.widthScale); // Assuming uniform width and length scale
                    newTree.transform.rotation = Quaternion.Euler(0, tree.rotation * 180f / Mathf.PI, 0); // Convert rotation to degrees

                    // Set to static
                    newTree.isStatic = true;

                    // Set rendering layers
                    newTree.layer = LayerMask.GetMask("Terrain");
                    if (newTree.TryGetComponent<MeshRenderer>(out MeshRenderer renderer))
                    {
                        renderer.renderingLayerMask |= config.renderingLayerMask;
                    }
                }
            }

            if (config.copyDetail)
            {
                Transform grassParent = new GameObject("Grass").transform;
                grassParent.parent = meshTerrain.transform;

                var terrainData = terrain.terrainData;
                var scaleX = terrainData.size.x / terrainData.detailWidth;
                var scaleZ = terrainData.size.z / terrainData.detailHeight;
                Debug.LogDebug("Detail Prototypes: " + terrainData.detailPrototypes.Length);


                for (int d = 0; d < terrainData.detailPrototypes.Length; d++)
                {
                    var detailPrototype = terrainData.detailPrototypes[d];
                    var detailLayer = terrainData.GetDetailLayer(0, 0, terrainData.detailWidth, terrainData.detailHeight, d);
                    float targetDensity = detailPrototype.density;
                    Debug.LogDebug("Target Coverage: " + detailPrototype.targetCoverage);
                    for (int x = 0; x < terrainData.detailWidth; x++)
                    {
                        for (int y = 0; y < terrainData.detailHeight; y++)
                        {
                            var layerDensity = detailLayer[y, x] / 255f;
                            float posX = x * scaleX + terrain.transform.position.x;
                            float posZ = y * scaleZ + terrain.transform.position.z;
                            float perlinNoise = Mathf.PerlinNoise(posX, posZ);
                            if (perlinNoise * layerDensity * targetDensity > 0.9f)
                            {
                                //Debug.Log("Density factor: " + perlinNoise * layerDensity * targetDensity);
                                var pos = new Vector3(posX, 0, posZ);
                                pos.y = terrain.SampleHeight(pos);
                                var detail = GameObject.Instantiate(terrainData.detailPrototypes[d].prototype, pos, Quaternion.Euler(0, UnityEngine.Random.Range(0, 359), 0), grassParent);

                                var scale = UnityEngine.Random.Range(detailPrototype.minWidth, detailPrototype.maxWidth);
                                var height = UnityEngine.Random.Range(detailPrototype.minHeight, detailPrototype.maxHeight);
                                detail.transform.localScale = new Vector3(scale, height, scale);
                            }

                        }
                    }
                }
            }

            return meshTerrain;
        }

        private static void GenerateMeshData(TerraMeshConfig config,
                                             MeshifyTerrainData terrainData)
        {
            // To avoid using terrainData methods, that are unavailable off the main thread
            float SampleHeightLocal(Vector3 pos)
            {
                float heightmapX = (pos.x - terrainData.terrainPosition.x) / terrainData.terrainStepX;
                float heightmapZ = (pos.z - terrainData.terrainPosition.z) / terrainData.terrainStepZ;
                int x = Mathf.Clamp(Mathf.RoundToInt(heightmapX), 0, terrainData.heightmapResolution - 1);
                int z = Mathf.Clamp(Mathf.RoundToInt(heightmapZ), 0, terrainData.heightmapResolution - 1);
                return terrainData.heightmapData[z, x] * terrainData.terrainHeight; // TODO check if order z,x is correct here
            }

            int actualCellStep;

            HashSet<Vector3> holeVertices = new HashSet<Vector3>();

            if (config.levelBounds != null) // Use the level bounds to determine the mesh density
            {
                terrainData.terrainBounds.center += terrainData.terrainPosition;
                float terrainSize = Mathf.Max(terrainData.terrainBounds.extents.x, terrainData.terrainBounds.extents.z);
                //Debug.LogDebug"Terrain center: " + terrainBounds.center + " Terrain Size: " + terrainSize);

                Vector3 levelCenter = config.levelBounds.Value.center;
                float levelSize = Mathf.Max(config.levelBounds.Value.extents.x, config.levelBounds.Value.extents.z);
                //Debug.LogDebug"Level Center: " + levelCenter + " Level Size: " + levelSize);

                int minMeshStep = config.minMeshStep;
                if (config.targetVertexCount > 0)
                {
                    minMeshStep = Mathf.CeilToInt(Mathf.Sqrt(levelSize * levelSize / (terrainData.terrainStepX * terrainData.terrainStepZ * config.targetVertexCount)));
                }

                //Debug.LogDebug"Base Density Factor: " + minMeshStep);

                QuadTree rootNode = new QuadTree(terrainData.terrainBounds);
                rootNode.Subdivide(config.levelBounds.Value, new Vector2(terrainData.terrainStepX, terrainData.terrainStepZ), minMeshStep,
                                    config.maxMeshStep, config.falloffSpeed, terrainSize - levelSize);

                // Generate vertices from 4 corners of leaf nodes of the quadtree
                HashSet<Vector3> uniqueVertices = new HashSet<Vector3>();
                List<QuadTree> leafNodes = new List<QuadTree>();
                rootNode.GetLeafNodes(leafNodes);

                foreach (var node in leafNodes)
                {
                    Vector3[] corners = new Vector3[]
                    {
                        new Vector3(node.bounds.min.x, 0, node.bounds.min.z),
                        new Vector3(node.bounds.min.x, 0, node.bounds.max.z),
                        new Vector3(node.bounds.max.x, 0, node.bounds.min.z),
                        new Vector3(node.bounds.max.x, 0, node.bounds.max.z)
                    };

                    foreach (var corner in corners)
                    {
                        if (uniqueVertices.Add(corner))
                        {
                            float height = SampleHeightLocal(corner);
                            Vector3 vertex = new Vector3(corner.x - terrainData.terrainPosition.x, height, corner.z - terrainData.terrainPosition.z);
                            terrainData.vertices.Add(vertex);

                            Vector2 uv = new Vector2(vertex.x / terrainData.terrainWidth, vertex.z / terrainData.terrainLength);
                            terrainData.uvs.Add(uv);

                            if (config.carveHoles)
                            {
                                int heightmapX = (int)((vertex.x) / terrainData.terrainStepX);
                                int heightmapZ = (int)((vertex.z) / terrainData.terrainStepZ);
                                heightmapX = Mathf.Clamp(heightmapX, 0, terrainData.holesResolution - 1);
                                heightmapZ = Mathf.Clamp(heightmapZ, 0, terrainData.holesResolution - 1);

                                // Check if the vertex is inside a terrain hole
                                if (!terrainData.holesData[heightmapX, heightmapZ])
                                {
                                    holeVertices.Add(vertex);
                                }
                            }
                        }
                    }
                }

                Debug.LogDebug("Sampled vertices: " + terrainData.vertices.Count);

                var polygon = new Polygon();

                for (int i = 0; i < terrainData.vertices.Count; i++)
                {
                    Vertex triNetVertex = terrainData.vertices[i].ToTriangleNetVertex(terrainData.uvs[i], 1);
                    polygon.Add(triNetVertex);
                }

                // Configure triangulation options
                int maxAdditionalVertices = Mathf.Min(terrainData.vertices.Count / 4, 30000);
                ConstraintOptions options = new ConstraintOptions() { ConformingDelaunay = false, SegmentSplitting = 2 };
                QualityOptions quality = new QualityOptions() { MinimumAngle = 20.0f, SteinerPoints = config.refineMesh ? maxAdditionalVertices : 0 };

                System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
                // Perform triangulation
                sw.Restart();
                var mesh2d = polygon.Triangulate(options, quality);
                sw.Stop();
                Debug.LogDebug("Triangulation time: " + sw.ElapsedMilliseconds + "ms");
                Debug.LogDebug("Final vertices: " + mesh2d.Vertices.Count);

                // Convert the 2D mesh to Unity mesh
                terrainData.vertices = new List<Vector3>();
                Dictionary<int, int> vertexIDs = new Dictionary<int, int>();
                terrainData.uvs = new List<Vector2>();

                foreach (ITriangle triangle in mesh2d.Triangles)
                {
                    int v0 = triangle.GetVertexID(0);
                    int v1 = triangle.GetVertexID(1);
                    int v2 = triangle.GetVertexID(2);

                    if (!vertexIDs.ContainsKey(v0))
                    {
                        vertexIDs[v0] = terrainData.vertices.Count;
                        terrainData.vertices.Add(triangle.GetVertex(0).ToVector3(1));
                        terrainData.uvs.Add(triangle.GetVertex(0).UV);
                    }
                    if (!vertexIDs.ContainsKey(v1))
                    {
                        vertexIDs[v1] = terrainData.vertices.Count;
                        terrainData.vertices.Add(triangle.GetVertex(1).ToVector3(1));
                        terrainData.uvs.Add(triangle.GetVertex(1).UV);
                    }
                    if (!vertexIDs.ContainsKey(v2))
                    {
                        vertexIDs[v2] = terrainData.vertices.Count;
                        terrainData.vertices.Add(triangle.GetVertex(2).ToVector3(1));
                        terrainData.uvs.Add(triangle.GetVertex(2).UV);
                    }

                    if (config.carveHoles)
                    {
                        if (holeVertices.Contains(terrainData.vertices[vertexIDs[v0]]) ||
                            holeVertices.Contains(terrainData.vertices[vertexIDs[v1]]) ||
                            holeVertices.Contains(terrainData.vertices[vertexIDs[v2]]))
                        {
                            continue;
                        }
                    }

                    terrainData.triangles.Add(vertexIDs[v0]);
                    terrainData.triangles.Add(vertexIDs[v2]);
                    terrainData.triangles.Add(vertexIDs[v1]);
                }
            }
            else // Uniform meshing if no level bounds are set
            {
                int minMeshStep = config.minMeshStep;
                // Calculate density factor to achieve target vertex count
                if (config.targetVertexCount > 0)
                {
                    minMeshStep = Mathf.CeilToInt(Mathf.Sqrt(terrainData.heightmapResolution * terrainData.heightmapResolution / config.targetVertexCount));
                }

                actualCellStep = Mathf.Max(minMeshStep, 1);
                //Debug.LogDebug"Density Factor: " + actualCellStep);

                // Calculate grid dimensions after applying density factor
                int gridWidth = Mathf.FloorToInt(terrainData.heightmapResolution / actualCellStep);
                int gridHeight = Mathf.FloorToInt(terrainData.heightmapResolution / actualCellStep);

                // Generate vertices
                for (int z = 0; z <= gridHeight; z++)
                {
                    for (int x = 0; x <= gridWidth; x++)
                    {
                        // Convert grid coordinates back to heightmap coordinates
                        int heightmapX = x * actualCellStep;
                        int heightmapZ = z * actualCellStep;

                        // Clamp to prevent accessing outside heightmap bounds
                        heightmapX = Mathf.Min(heightmapX, terrainData.heightmapResolution - 1);
                        heightmapZ = Mathf.Min(heightmapZ, terrainData.heightmapResolution - 1);

                        float height = terrainData.heightmapData[heightmapZ, heightmapX] * terrainData.terrainHeight; // TODO check if order z,x is correct here
                        Vector3 vertex = new Vector3(heightmapX * terrainData.terrainStepX, height, heightmapZ * terrainData.terrainStepZ);
                        terrainData.vertices.Add(vertex);

                        Vector2 uv = new Vector2(heightmapX * terrainData.uvStepX, heightmapZ * terrainData.uvStepZ);
                        terrainData.uvs.Add(uv);

                        if (config.carveHoles)
                        {
                            heightmapX = Mathf.Clamp(heightmapX, 0, terrainData.holesResolution - 1);
                            heightmapZ = Mathf.Clamp(heightmapZ, 0, terrainData.holesResolution - 1);

                            // Check if the vertex is inside a terrain hole
                            if (!terrainData.holesData[heightmapX, heightmapZ])
                            {
                                holeVertices.Add(vertex);
                            }
                        }
                    }
                }

                Debug.LogDebug("Sampled vertices: " + terrainData.vertices.Count);

                // Generate triangles using grid coordinates
                for (int z = 0; z < gridHeight; z++)
                {
                    for (int x = 0; x < gridWidth; x++)
                    {
                        // Calculate vertex indices in the grid
                        int vertexIndex = z * (gridWidth + 1) + x;

                        if (config.carveHoles)
                        {
                            if (holeVertices.Contains(terrainData.vertices[vertexIndex]) ||
                                holeVertices.Contains(terrainData.vertices[vertexIndex + 1]) ||
                                holeVertices.Contains(terrainData.vertices[vertexIndex + (gridWidth + 1)]) ||
                                holeVertices.Contains(terrainData.vertices[vertexIndex + (gridWidth + 1) + 1]))
                            {
                                continue;
                            }
                        }

                        // First triangle
                        terrainData.triangles.Add(vertexIndex);                     // Current vertex
                        terrainData.triangles.Add(vertexIndex + (gridWidth + 1));   // Vertex below
                        terrainData.triangles.Add(vertexIndex + (gridWidth + 1) + 1); // Vertex below and right

                        // Second triangle
                        terrainData.triangles.Add(vertexIndex);                     // Current vertex
                        terrainData.triangles.Add(vertexIndex + (gridWidth + 1) + 1); // Vertex below and right
                        terrainData.triangles.Add(vertexIndex + 1);                 // Vertex to the right
                    }
                }
            }
        }

#if UNITY_EDITOR
        public static Texture2D[] GetSplatmapsAsTextures(this Terrain terrain, string? savePath = null)
#else
        public static Texture2D[] GetSplatmapsAsTextures(this Terrain terrain)
#endif
        {

            TerrainData terrainData = terrain.terrainData;
            int alphamapWidth = terrainData.alphamapWidth;
            int alphamapHeight = terrainData.alphamapHeight;
            float[,,] splatmapData = terrainData.GetAlphamaps(0, 0, alphamapWidth, alphamapHeight);
            int numSplatmaps = terrainData.alphamapLayers;

            // Calculate how many RGBA textures we need
            int numTextures = Mathf.CeilToInt(numSplatmaps / 4f);
            Texture2D[] splatmaps = new Texture2D[numTextures];

            for (int textureIndex = 0; textureIndex < numTextures; textureIndex++)
            {
                Texture2D packedTexture = new Texture2D(alphamapWidth, alphamapHeight, TextureFormat.RGBA32, false);
                Color[] packedColors = new Color[alphamapWidth * alphamapHeight];

                for (int y = 0; y < alphamapHeight; y++)
                {
                    for (int x = 0; x < alphamapWidth; x++)
                    {
                        float r = textureIndex * 4 + 0 < numSplatmaps ? splatmapData[y, x, textureIndex * 4 + 0] : 0f;
                        float g = textureIndex * 4 + 1 < numSplatmaps ? splatmapData[y, x, textureIndex * 4 + 1] : 0f;
                        float b = textureIndex * 4 + 2 < numSplatmaps ? splatmapData[y, x, textureIndex * 4 + 2] : 0f;
                        float a = textureIndex * 4 + 3 < numSplatmaps ? splatmapData[y, x, textureIndex * 4 + 3] : 0f;

                        packedColors[y * alphamapWidth + x] = new Color(r, g, b, a);
                    }
                }

                packedTexture.SetPixels(packedColors);
                packedTexture.Apply();

#if UNITY_EDITOR
                if (savePath != null)
                {
                    string filePath = Path.Combine(savePath, $"Splatmap_{textureIndex}_{terrain.name}.png");
                    byte[] bytes = packedTexture.EncodeToPNG();

                    File.WriteAllBytes(filePath, bytes);  // Write the texture data to disk

                    AssetDatabase.ImportAsset(filePath); // Import the asset into the project

                    // Retrieve the texture asset directly
                    splatmaps[textureIndex] = AssetDatabase.LoadAssetAtPath<Texture2D>(filePath);
                }
                else
                {
                    splatmaps[textureIndex] = packedTexture;
                }
#else
                splatmaps[textureIndex] = packedTexture;
#endif

                Debug.LogDebug($"Saved splatmap_{textureIndex}.png with layers:" +
                    $"\nR: Layer {textureIndex * 4 + 0}" +
                    (textureIndex * 4 + 1 < numSplatmaps ? $"\nG: Layer {textureIndex * 4 + 1}" : "") +
                    (textureIndex * 4 + 2 < numSplatmaps ? $"\nB: Layer {textureIndex * 4 + 2}" : "") +
                    (textureIndex * 4 + 3 < numSplatmaps ? $"\nA: Layer {textureIndex * 4 + 3}" : ""));
            }

            return splatmaps;
        }


        /// <summary>
        /// Sets up a material with terrain textures based on the Terrain object.
        /// </summary>
        /// <param name="targetMaterial">The Material to set up.</param>
        /// <param name="terrain">The Terrain object to get textures from.</param>
        /// <param name="savePath">The path to save the splatmaps as PNG files. If null, the splatmaps will not be saved.</param>
        /// <remarks>
        /// This method sets up the material with textures from the Terrain object's TerrainLayers.
        /// The material should use MeshTerrainLit shader.
        /// The splatmaps are saved as PNG files if the savePath is provided.
        /// </remarks>
#if UNITY_EDITOR
        public static void SetupMaterialFromTerrain(this Material targetMaterial, Terrain terrain, string? savePath = null)
#else
        public static void SetupMaterialFromTerrain(this Material targetMaterial, Terrain terrain)
#endif
        {
            TerrainLayer[] terrainLayers = terrain.terrainData.terrainLayers;
            int layerCount = terrainLayers.Length;

            if (layerCount > 8)
            {
                Debug.LogWarning("Terrain has more than 8 layers. Only the first 8 will be used.");
                layerCount = 8;
            }

            // Get splatmaps
#if UNITY_EDITOR
            Texture2D[] splatmaps = terrain.GetSplatmapsAsTextures(savePath);
#else
            Texture2D[] splatmaps = terrain.GetSplatmapsAsTextures();
#endif

            for (int i = 0; i < splatmaps.Length; i++)
            {
                targetMaterial.SetTexture($"_Splatmap_{i}", splatmaps[i]);
            }

            // Process each layer
            for (int i = 0; i < layerCount; i++)
            {
                TerrainLayer layer = terrainLayers[i];

                //Skip empty layers
                if (layer.diffuseTexture == null)
                {
                    continue;
                }

                // Set textures
                targetMaterial.SetTexture($"_Albedo_{i}", layer.diffuseTexture);
                targetMaterial.SetTexture($"_Normals_{i}", layer.normalMapTexture);
                targetMaterial.SetTexture($"_Mask_{i}", layer.maskMapTexture);

                targetMaterial.SetColor($"_Color_Tint_{i}", layer.diffuseRemapMax);
                targetMaterial.SetFloat($"_Normal_Scale_{i}", layer.normalScale);

                // Set tiling and offset
                targetMaterial.SetVector($"_Tiling_{i}", new Vector2(layer.tileSize.x, layer.tileSize.y));
                targetMaterial.SetVector($"_Offset_{i}", new Vector2(layer.tileOffset.x, layer.tileOffset.y));

                // Set remapping values
                if (layer.maskMapTexture == null)
                {
                    targetMaterial.SetVector($"_Metallic_Remapping_{i}", new Vector2(layer.metallic, layer.metallic));
                    targetMaterial.SetVector($"_AO_Remapping_{i}", new Vector2(1, 1));
                    if (GraphicsFormatUtility.GetAlphaComponentCount(layer.diffuseTexture.format) > 0)
                    {
                        targetMaterial.SetVector($"_Smoothness_Remapping_{i}", new Vector2(0, 1));
                    }
                    else
                    {
                        targetMaterial.SetVector($"_Smoothness_Remapping_{i}", new Vector2(layer.smoothness, layer.smoothness));
                    }
                }
                else
                {
                    targetMaterial.SetVector($"_Metallic_Remapping_{i}", new Vector2(layer.maskMapRemapMin.x, layer.maskMapRemapMax.x));
                    targetMaterial.SetVector($"_AO_Remapping_{i}", new Vector2(layer.maskMapRemapMin.y, layer.maskMapRemapMax.y));
                    targetMaterial.SetVector($"_Smoothness_Remapping_{i}", new Vector2(layer.maskMapRemapMin.w, layer.maskMapRemapMax.w));
                }
            }

            // Clear unused layers
            for (int i = layerCount; i < 8; i++)
            {
                Texture2D emptyTex = Texture2D.blackTexture;

                targetMaterial.SetTexture($"_Albedo_{i}", emptyTex);
                targetMaterial.SetTexture($"_Normals_{i}", emptyTex);
                targetMaterial.SetTexture($"_Mask_{i}", emptyTex);
                targetMaterial.SetColor($"_Color_Tint_{i}", Vector4.zero);
                targetMaterial.SetFloat($"_Normal_Scale_{i}", 0f);
                targetMaterial.SetVector($"_Tiling_{i}", Vector2.one);
                targetMaterial.SetVector($"_Offset_{i}", Vector2.zero);
                targetMaterial.SetVector($"_Metallic_Remapping_{i}", new Vector2(0, 1));
                targetMaterial.SetVector($"_AO_Remapping_{i}", new Vector2(0, 1));
                targetMaterial.SetVector($"_Smoothness_Remapping_{i}", new Vector2(0, 1));
            }
        }
    }

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
    public class MeshifyTerrainData
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

        public List<Vector3> vertices;
        public List<int> triangles;
        public List<Vector2> uvs;

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