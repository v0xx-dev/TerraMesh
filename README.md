# TerraMesh
This is a simple plugin that converts Unity Terrain to mesh objects. Based on a fork of Triangle.NET (same license applies to those parts of the code).

## How to use:

### In Unity Editor

- Select "Tools/TerraMesh" to show the main window, hover over the variable names for more information.

### At runtime

1. Create an instance of meshing configuration TerraMeshConfig
2. Use .Meshify method on the desired Unity Terrain object to turn it into mesh
3. Obtain the material instance from the mesh renderer of newly created mesh terrain and run .SetupMaterialFromTerrain method to copy texture references from TerrainLayers

More information here: https://discord.com/channels/1168655651455639582/1303914349533990983 

### Information

It has two ways of meshing:
1. **Uniform**: if you don't provide a box collider it will apply a regular grid with a specified step size (1 = exactly the same as original, 2 = resolution halved, 3 = resolution/3, etc).
2. **Variable**: Provide a box collider intersecting with your terrain and apply non-uniform meshing: within the box it will use min step size, outside the box the step will gradually increase towards the edges (only in multiples of two).
*HIGHLY recommend to use the 2nd option!*

- **Falloff factor** regulates how fast the step size changes.
- Use **refine mesh** option if you see very thin triangles appearing, this can happen if your terrain is rectangular or a very high falloff used with big maximum step, realistically this should be disabled in most cases.
- **Copy trees** toogle will save trees from your terrain as separate objects
- **Carve holes** toggle will carve holes in the mesh if they exist in your terrain. 

Grass is NOT supported, probably won't be.

## **Performance considerations!**

- Set your mesh terrain/trees objects to **static** and set the material of your trees/grass to use "GPU Instanced" option - this will allow batching, instead of rendering each object in separate draw calls.
 - To keep rendering performance the same, your terrain should have less than 100k vertices, for vanilla sized moons the number is close to **20-30k**.
- It automatically applies a mesh collider so be careful with high polygon counts! **You can run the tool with lower quality and use the simpler mesh for the collider!**
- Ideally, have a high res terrain for your playable area and low res outer terrain (w/o collider) to create a feeling of a vast environment.