This is a simple plugin that converts Unity Terrain to mesh objects. Based on a fork of Triangle.NET (same license applies to those parts of the code).

# How to use:

## In Unity Editor

- Select "Tools/TerraMesh" to show the main window, hover over the variable names for more information.

## At runtime

1. Create an instance of meshing configuration TerraMeshConfig
2. Use .Meshify method on the desired Unity Terrain object to turn it into mesh
3. Obtain the material instance from the mesh renderer of newly created mesh terrain and run .SetupMaterialFromTerrain method to copy texture references from TerrainLayers

More information here: https://discord.com/channels/1168655651455639582/1303914349533990983 