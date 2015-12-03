using System;
using System.IO;
using System.Reflection;
using Illumina.Common;

/// <summary>
///     Contains extension methods for the StreamWriter class
/// </summary>
public static class StreamWriterExtension
{
    /// <summary>
    ///     Gets the current read position of the StreamWriter.
    /// </summary>
    public static long TruePosition(this StreamWriter streamWriter)
    {
        Type writerType = streamWriter.GetType();

        int bufferPos = 0;

        if (CrossPlatform.IsThisMono())
        {
            // retrieve the character position
            bufferPos = (int) writerType.InvokeMember(
                "charPos",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamWriter,
                null);
        }
        else
        {
            // retrieve the character position
            bufferPos = (int) writerType.InvokeMember(
                "charPos",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamWriter,
                null);
        }

        // using both we can calculate the true file position
        return streamWriter.BaseStream.Position + bufferPos;
    }
}