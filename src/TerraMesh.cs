using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using System.IO;
using UnityEngine.Experimental.Rendering;

namespace TriangleNet.Unity
{
    public class TerraMesh : EditorWindow
{
    private Terrain terrain;
    private bool refineMesh = false;
    private bool carveHoles = false;
    private bool copyTrees = false;
    private bool copyDetail = false;
    
    private BoxCollider levelBounds;
    private int targetVertexCount = -1;
    private int minMeshStep = 1;
    private int maxMeshStep = 16;
    private float falloffSpeed = 2f;

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
                string currentTime = DateTime.Now.ToString("ddMHHmmss"); 
                string terrainFolderPath = Path.Combine(saveFolderPath, terrain.name + "_" + currentTime); // Subfolder for the terrain

                // Create subfolder if it doesn't exist
                Directory.CreateDirectory(terrainFolderPath);

                Material terrainMaterial = new Material(Shader.Find("Shader Graphs/MeshTerrainLit"));
                string materialPath = Path.Combine(terrainFolderPath, terrain.name + "_Lit.mat");
                AssetDatabase.CreateAsset(terrainMaterial, materialPath);

                meshTerrain = terrain.Meshify(terrainMaterial, levelBounds, targetVertexCount, minMeshStep,
                                                maxMeshStep, falloffSpeed, refineMesh, carveHoles, copyTrees, copyDetail);

                terrain.SetupMaterialFromTerrain(terrainMaterial, true, terrainFolderPath);

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

    static class TerraMeshExtensions
    {
        
        public static GameObject Meshify(this Terrain terrain, Material terrainMaterial, BoxCollider levelBounds,
                                        int targetVertexCount, float minCellStep, float maxCellStep, float falloffSpeed,
                                        bool refineMesh = false, bool carveHoles = false, bool copyTrees = false, bool copyDetail = false)
        {
            GameObject meshTerrain = new GameObject("MeshTerrain_" + terrain.name);
            MeshFilter meshFilter = meshTerrain.AddComponent<MeshFilter>();
            MeshRenderer meshRenderer = meshTerrain.AddComponent<MeshRenderer>();
            MeshCollider meshCollider = meshTerrain.AddComponent<MeshCollider>();

            Mesh mesh = new Mesh();
            mesh.name = "MeshTerrain_" + terrain.name;

            List<Vector3> vertices = new List<Vector3>();
            List<Vector2> uvs = new List<Vector2>();
            List<int> triangles = new List<int>();

            meshTerrain.transform.position = terrain.transform.position;

            float terrainWidth = terrain.terrainData.size.x;
            float terrainLength = terrain.terrainData.size.z;
            float terrainHeight = terrain.terrainData.size.y;

            float terrainStepX = terrainWidth / (terrain.terrainData.heightmapResolution - 1);
            float terrainStepZ = terrainLength / (terrain.terrainData.heightmapResolution - 1);
            
            float uvStepX = 1.0f / (terrain.terrainData.heightmapResolution - 1);
            float uvStepZ = 1.0f / (terrain.terrainData.heightmapResolution - 1);
            
            int actualCellStep;

            HashSet<Vector3> holeVertices = new HashSet<Vector3>();

            if (levelBounds != null)
            {
                Bounds terrainBounds = terrain.terrainData.bounds;
                terrainBounds.center += terrain.transform.position;
                float terrainSize = Mathf.Max(terrainBounds.extents.x, terrainBounds.extents.z);
                //Debug.Log("Terrain center: " + terrainBounds.center + " Terrain Size: " + terrainSize);

                Vector3 levelCenter = levelBounds.bounds.center;
                float levelSize = Mathf.Max(levelBounds.bounds.extents.x, levelBounds.bounds.extents.z);
                //Debug.Log("Level Center: " + levelCenter + " Level Size: " + levelSize);

                if (targetVertexCount > 0)
                {
                    minCellStep = Mathf.Sqrt((levelSize * levelSize) / (terrainStepX * terrainStepZ * targetVertexCount));
                    minCellStep = Mathf.CeilToInt(minCellStep);
                }

                //Debug.Log("Base Density Factor: " + minCellStep);

                QuadTree rootNode = new QuadTree(terrainBounds);
                rootNode.Subdivide(levelBounds.bounds, new Vector2(terrainStepX, terrainStepZ), minCellStep,
                                    maxCellStep, falloffSpeed, terrainSize - levelSize);

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
                            float height = terrain.SampleHeight(corner);
                            Vector3 vertex = new Vector3(corner.x - terrain.transform.position.x, height, corner.z - terrain.transform.position.z);
                            vertices.Add(vertex);

                            Vector2 uv = new Vector2(vertex.x / terrainWidth, vertex.z / terrainLength);
                            uvs.Add(uv);

                            if (carveHoles)
                            {
                                int heightmapX = (int)((vertex.x) / terrainStepX );
                                int heightmapZ = (int)((vertex.z) / terrainStepZ);
                                heightmapX = Mathf.Clamp(heightmapX, 0, terrain.terrainData.holesResolution - 1);
                                heightmapZ = Mathf.Clamp(heightmapZ, 0, terrain.terrainData.holesResolution - 1);

                                // Check if the vertex is inside a terrain hole
                                if (!terrain.terrainData.GetHoles(heightmapX, heightmapZ, 1, 1)[0, 0])
                                {
                                    holeVertices.Add(vertex);
                                }
                            }
                        }
                    }
                }

                Debug.Log("Sampled vertices: " + vertices.Count);

                var polygon = new Polygon();

                for (int i = 0; i < vertices.Count; i++)
                {
                    Vertex triNetVertex = vertices[i].ToTriangleNetVertex(uvs[i], 1);
                    polygon.Add(triNetVertex);
                }

                // Configure triangulation options
                ConstraintOptions options = new ConstraintOptions() { ConformingDelaunay = false, SegmentSplitting = 2};
                QualityOptions quality = new QualityOptions() { MinimumAngle = 20.0f, SteinerPoints = refineMesh ? -1 : 0};

                System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
                // Perform triangulation
                sw.Restart();
                var mesh2d = polygon.Triangulate(options, quality);
                sw.Stop();
                Debug.Log("Triangulation time: " + sw.ElapsedMilliseconds + "ms");
                Debug.Log("Final vertices: " + mesh2d.Vertices.Count);

                // Convert the 2D mesh to Unity mesh
                vertices = new List<Vector3>();
                Dictionary<int, int> vertexIDs = new Dictionary<int, int>();
                uvs = new List<Vector2>();

                foreach (ITriangle triangle in mesh2d.Triangles)
                {
                    int v0 = triangle.GetVertexID(0);
                    int v1 = triangle.GetVertexID(1);
                    int v2 = triangle.GetVertexID(2);

                    if (!vertexIDs.ContainsKey(v0))
                    {
                        vertexIDs[v0] = vertices.Count;
                        vertices.Add(triangle.GetVertex(0).ToVector3(1));
                        uvs.Add(triangle.GetVertex(0).UV);
                    }
                    if (!vertexIDs.ContainsKey(v1))
                    {
                        vertexIDs[v1] = vertices.Count;
                        vertices.Add(triangle.GetVertex(1).ToVector3(1));
                        uvs.Add(triangle.GetVertex(1).UV);
                    }
                    if (!vertexIDs.ContainsKey(v2))
                    {
                        vertexIDs[v2] = vertices.Count;
                        vertices.Add(triangle.GetVertex(2).ToVector3(1));
                        uvs.Add(triangle.GetVertex(2).UV);
                    }

                    if (carveHoles)
                    {
                        if (holeVertices.Contains(vertices[vertexIDs[v0]]) ||
                            holeVertices.Contains(vertices[vertexIDs[v1]]) ||
                            holeVertices.Contains(vertices[vertexIDs[v2]]))
                        {
                            continue;
                        }
                    }

                    triangles.Add(vertexIDs[v0]);
                    triangles.Add(vertexIDs[v2]);
                    triangles.Add(vertexIDs[v1]);
                }
            }
            else // Uniform meshing if no level bounds are set
            {
                // Calculate density factor to achieve target vertex count
                if (targetVertexCount > 0)
                {
                    minCellStep = Mathf.Sqrt(terrain.terrainData.heightmapResolution * terrain.terrainData.heightmapResolution / targetVertexCount);
                    minCellStep = Mathf.CeilToInt(minCellStep);
                }

                actualCellStep = (int)Mathf.Max(minCellStep, 1);
                //Debug.Log("Density Factor: " + actualCellStep);

                // Calculate grid dimensions after applying density factor
                int gridWidth = Mathf.FloorToInt((terrain.terrainData.heightmapResolution) / actualCellStep);
                int gridHeight = Mathf.FloorToInt((terrain.terrainData.heightmapResolution) / actualCellStep);

                // Generate vertices
                for (int z = 0; z <= gridHeight; z++)
                {
                    for (int x = 0; x <= gridWidth; x++)
                    {
                        // Convert grid coordinates back to heightmap coordinates
                        int heightmapX = x * actualCellStep;
                        int heightmapZ = z * actualCellStep;

                        // Clamp to prevent accessing outside heightmap bounds
                        heightmapX = Mathf.Min(heightmapX, terrain.terrainData.heightmapResolution - 1);
                        heightmapZ = Mathf.Min(heightmapZ, terrain.terrainData.heightmapResolution - 1);

                        float height = terrain.terrainData.GetHeight(heightmapX, heightmapZ);
                        Vector3 vertex = new Vector3(heightmapX * terrainStepX, height, heightmapZ * terrainStepZ);
                        vertices.Add(vertex);

                        Vector2 uv = new Vector2(heightmapX * uvStepX, heightmapZ * uvStepZ);
                        uvs.Add(uv);

                        if (carveHoles)
                        {
                            heightmapX = Mathf.Clamp(heightmapX, 0, terrain.terrainData.holesResolution - 1);
                            heightmapZ = Mathf.Clamp(heightmapZ, 0, terrain.terrainData.holesResolution - 1);

                            // Check if the vertex is inside a terrain hole
                            if (!terrain.terrainData.GetHoles(heightmapX, heightmapZ, 1, 1)[0, 0])
                            {
                                holeVertices.Add(vertex);
                            }
                        }
                    }
                }

                Debug.Log("Sampled vertices: " + vertices.Count);

                // Generate triangles using grid coordinates
                for (int z = 0; z < gridHeight; z++)
                {
                    for (int x = 0; x < gridWidth; x++)
                    {
                        // Calculate vertex indices in the grid
                        int vertexIndex = z * (gridWidth + 1) + x;

                        if (carveHoles)
                        {
                            if (holeVertices.Contains(vertices[vertexIndex]) ||
                                holeVertices.Contains(vertices[vertexIndex + 1]) ||
                                holeVertices.Contains(vertices[vertexIndex + (gridWidth + 1)]) ||
                                holeVertices.Contains(vertices[vertexIndex + (gridWidth + 1) + 1]))
                            {
                                continue;
                            }
                        }

                        // First triangle
                        triangles.Add(vertexIndex);                     // Current vertex
                        triangles.Add(vertexIndex + (gridWidth + 1));   // Vertex below
                        triangles.Add(vertexIndex + (gridWidth + 1) + 1); // Vertex below and right

                        // Second triangle
                        triangles.Add(vertexIndex);                     // Current vertex
                        triangles.Add(vertexIndex + (gridWidth + 1) + 1); // Vertex below and right
                        triangles.Add(vertexIndex + 1);                 // Vertex to the right
                    }
                }
            }

            mesh.indexFormat = vertices.Count > 65535 ? UnityEngine.Rendering.IndexFormat.UInt32 : UnityEngine.Rendering.IndexFormat.UInt16;
            mesh.vertices = vertices.ToArray();
            mesh.uv = uvs.ToArray();
            mesh.triangles = triangles.ToArray();
            mesh.Optimize(); // Clean unused vertices
            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            mesh.RecalculateTangents();

            meshFilter.mesh = mesh;
            meshCollider.sharedMesh = mesh;
            meshRenderer.material = terrainMaterial;

            if (copyTrees)
            {
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

                    // Calculate tree position
                    Vector3 treeWorldPos = Vector3.Scale(tree.position, terrain.terrainData.size) + terrain.transform.position;

                    // Set tree transform
                    newTree.transform.position = treeWorldPos;
                    newTree.transform.localScale = new Vector3(tree.widthScale, tree.heightScale, tree.widthScale); // Assuming uniform width and length scale
                    newTree.transform.rotation = Quaternion.Euler(0, tree.rotation * 180f/Mathf.PI, 0); // Convert rotation to degrees
                }
            }

            if (copyDetail)
            {
                Transform grassParent = new GameObject("Grass").transform;
                grassParent.parent = meshTerrain.transform;

                var terrainData = terrain.terrainData;
                var scaleX = terrainData.size.x / terrainData.detailWidth;
                var scaleZ = terrainData.size.z / terrainData.detailHeight;
                Debug.Log("Detail Prototypes: " + terrainData.detailPrototypes.Length);
                
                
                for (int d = 0; d < terrainData.detailPrototypes.Length; d++)
                {
                    var detailPrototype = terrainData.detailPrototypes[d];
                    var detailLayer = terrainData.GetDetailLayer(0, 0, terrainData.detailWidth, terrainData.detailHeight, d);
                    float targetDensity = detailPrototype.density;
                    Debug.Log("Target Coverage: " + detailPrototype.targetCoverage);
                    for (int x = 0; x < terrainData.detailWidth; x++)
                    {
                        for (int y = 0; y < terrainData.detailHeight; y++)
                        {
                            var layerDensity = detailLayer[y, x]/255f;
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

        public static Texture2D[] GetSplatmapsAsTextures(this Terrain terrain, string folderPath = "Assets/", bool save = false)
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

                if (save)
                {
                    string filePath = Path.Combine(folderPath, $"Splatmap_{textureIndex}_{terrain.name}.png");
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

                Debug.Log($"Saved splatmap_{textureIndex}.png with layers:" +
                    $"\nR: Layer {textureIndex * 4 + 0}" +
                    (textureIndex * 4 + 1 < numSplatmaps ? $"\nG: Layer {textureIndex * 4 + 1}" : "") +
                    (textureIndex * 4 + 2 < numSplatmaps ? $"\nB: Layer {textureIndex * 4 + 2}" : "") +
                    (textureIndex * 4 + 3 < numSplatmaps ? $"\nA: Layer {textureIndex * 4 + 3}" : ""));
            }

            return splatmaps;
        }

        public static void SetupMaterialFromTerrain(this Terrain terrain, Material targetMaterial, bool saveSplatmaps = false, string folderPath = "Assets/")
        {
            TerrainLayer[] terrainLayers = terrain.terrainData.terrainLayers;
            int layerCount = terrainLayers.Length;

            if (layerCount > 8)
            {
                Debug.LogWarning("Terrain has more than 8 layers. Only the first 8 will be used.");
                layerCount = 8;
            }

            // Get splatmaps
            Texture2D[] splatmaps = terrain.GetSplatmapsAsTextures(folderPath, saveSplatmaps);

            for (int i = 0; i < splatmaps.Length; i++)
            {
                targetMaterial.SetTexture($"_Splatmap_{i}", splatmaps[i]);
            }

            // Process each layer
            for (int i = 0; i < layerCount; i++)
            {
                TerrainLayer layer = terrainLayers[i];
                
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
    
    internal class QuadTree
    {
        public Bounds bounds;
        public QuadTree[] children;
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
                foreach (var child in children)
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
                foreach (var child in children)
                {
                    child.GetLeafNodes(leafNodes);
                }
            }
        } 
    }
}