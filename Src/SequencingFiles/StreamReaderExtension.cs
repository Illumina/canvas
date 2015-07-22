using System;
using System.IO;
using System.Reflection;
using Illumina.Common;

/// <summary>
///     Contains extension methods for the StreamReader class
/// </summary>
public static class StreamReaderExtension
{
    /// <summary>
    ///     Gets the current read position of the StreamReader.
    /// </summary>
    public static long TruePosition(this StreamReader streamReader)
    {
        Type readerType = streamReader.GetType();

        int bufferSize = 0;
        int bufferPos = 0;

        if (CrossPlatform.IsThisMono())
        {
            // retrieve the character position
            bufferPos = (int)readerType.InvokeMember(
                "pos",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamReader,
                null);

            // retrieve the buffer size
            bufferSize = (int)readerType.InvokeMember(
                "decoded_count",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamReader,
                null);
        }
        else
        {
            // retrieve the character position
            bufferPos = (int)readerType.InvokeMember(
                "charPos",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamReader,
                null);

            // retrieve the buffer size
            bufferSize = (int)readerType.InvokeMember(
                "charLen",
                BindingFlags.DeclaredOnly | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance |
                BindingFlags.GetField,
                null,
                streamReader,
                null);
        }

        // using both we can calculate the true file position
        return streamReader.BaseStream.Position - bufferSize + bufferPos;
    }

    /// <summary>
    ///     Sets the current read position of the StreamReader.
    /// </summary>
    public static void TrueSeek(this StreamReader streamReader, long offset, SeekOrigin origin)
    {
        streamReader.BaseStream.Seek(offset, origin);
        streamReader.DiscardBufferedData();
    }
}