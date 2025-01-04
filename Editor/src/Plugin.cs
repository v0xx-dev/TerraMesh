using System.Reflection;
using BepInEx;
using BepInEx.Logging;
using UnityEngine;
using System.IO;

namespace TerraMesh
{
    [BepInPlugin(PluginInfo.PLUGIN_GUID, PluginInfo.PLUGIN_NAME, PluginInfo.PLUGIN_VERSION)]
    public class TerraMeshPlugin : BaseUnityPlugin
    {
        public static TerraMeshPlugin Instance;
        internal static ManualLogSource StaticLogger;
        internal static Shader? terraMeshShader;

        private void Awake()
        {
                Instance = this;
                StaticLogger = Logger; 
                // Load the shader from an asset bundle
                string dllPath = Assembly.GetExecutingAssembly().Location;
                string bundlePath = Path.Combine(Path.GetDirectoryName(dllPath), "terramesh.assetbundle");
                var assetBundle = AssetBundle.LoadFromFile(bundlePath);
                terraMeshShader = assetBundle.LoadAsset<Shader>("MeshTerrainLit");
                if (terraMeshShader == null)
                {
                    Debug.LogError("Failed to load TerraMeshShader from asset bundle");
                }
                Logger.LogInfo($"Plugin {PluginInfo.PLUGIN_GUID} is loaded!");

        #if DEBUG
                // disable overhead of stack trace in dev build
                Application.SetStackTraceLogType(LogType.Log, StackTraceLogType.None);
                Application.SetStackTraceLogType(LogType.Warning, StackTraceLogType.None);
                Application.SetStackTraceLogType(LogType.Error, StackTraceLogType.None);
                Application.SetStackTraceLogType(LogType.Assert, StackTraceLogType.None);
        #endif
        }
    }

    public static class Debug
    {
        private static ManualLogSource Logger => TerraMeshPlugin.StaticLogger;
#if UNITY_EDITOR
        public static void Log(string message) => UnityEngine.Debug.Log(message);
        public static void LogError(string message) => UnityEngine.Debug.LogError(message);
        public static void LogWarning(string message) => UnityEngine.Debug.LogWarning(message);
        public static void LogDebug(string message) => UnityEngine.Debug.Log(message);
        public static void LogMessage(string message) => UnityEngine.Debug.Log(message);
        public static void LogFatal(string message) => UnityEngine.Debug.LogError(message);
#else
        public static void Log(string message) => Logger.LogInfo(message);
        public static void LogError(string message) => Logger.LogError(message);
        public static void LogWarning(string message) => Logger.LogWarning(message);
        public static void LogDebug(string message) => Logger.LogDebug(message);
        public static void LogMessage(string message) => Logger.LogMessage(message);
        public static void LogFatal(string message) => Logger.LogFatal(message);
#endif
    }
}