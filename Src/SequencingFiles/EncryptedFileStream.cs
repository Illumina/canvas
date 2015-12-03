using System;
using System.IO;
using System.Security.Cryptography;
using System.Text;

namespace SequencingFiles
{
    /// <summary>
    ///     The EncryptedFileStream uses a MemoryStream to hold the unencrypted contents.
    ///     When opening an encrypted file for reading, the file is decrypted and kept in
    ///     the memory stream. When writing to a file, the MemoryStream is encrypted when
    ///     the stream is closed.
    /// </summary>
    public class EncryptedFileStream : Stream
    {
        #region members

        private const string Salt = "IlluminaIsis";
        private readonly MemoryStream _baseMemoryStream;
        private readonly string _filePath;
        private readonly string _key;
        private FileStream _baseFileStream;
        private bool _isDisposed;
        private bool _isWriting;

        #endregion

        // constructor
        public EncryptedFileStream(string filePath, string key, FileAccess fileAccess)
        {
            _filePath = filePath;
            _key = key;
            _baseMemoryStream = new MemoryStream();

            Open(filePath, fileAccess);
        }

        #region Access properties

        /// <summary>
        ///     Returns true of this stream can be read from, false otherwise
        /// </summary>
        public override bool CanRead
        {
            get { return !_isWriting; }
        }

        /// <summary>
        ///     Returns false.
        /// </summary>
        public override bool CanSeek
        {
            get { return false; }
        }

        /// <summary>
        ///     Returns true if this stream is writable, false otherwise
        /// </summary>
        public override bool CanWrite
        {
            get { return _isWriting; }
        }

        #endregion

        #region Destructor & IDispose stuff

        /// <summary>
        ///     Destroys this instance
        /// </summary>
        ~EncryptedFileStream()
        {
            Dispose(false);
        }

        /// <summary>
        ///     Releases managed and optionally unmanaged resources
        /// </summary>
        /// <param name="disposing">True releases both managed and unmanaged, false releases only unmanaged resources</param>
        protected override void Dispose(bool disposing)
        {
            if (!_isDisposed)
            {
                if (disposing)
                {
                    if (_isWriting) Encrypt();
                    _baseMemoryStream.Close();
                    _baseFileStream.Close();
                }
                _isDisposed = true;
            }
            base.Dispose(disposing);
        }

        #endregion

        #region Properties

        /// <summary>
        ///     Gets the length in bytes of the stream.
        /// </summary>
        public override long Length
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        ///     Gets or sets the current position of this stream.
        /// </summary>
        public override long Position
        {
            get { throw new NotImplementedException(); }
            set { throw new NotImplementedException(); }
        }

        #endregion

        #region Methods

        /// <summary>
        ///     Returns true if the two hashes are equal
        /// </summary>
        private bool AreHashesEqual(byte[] hash1, byte[] hash2)
        {
            // sanity check: detect null arrays
            if ((hash1 == null) || (hash2 == null)) return false;

            // sanity check: check lengths
            if (hash1.Length != hash2.Length) return false;

            for (int i = 0; i < hash1.Length; ++i)
            {
                if (hash1[i] != hash2[i]) return false;
            }

            return true;
        }

        /// <summary>
        ///     fills the memory stream with the decrypted contents from the file stream
        /// </summary>
        private void Decrypt()
        {
            // read the hash from the file stream
            SHA256 shaM = new SHA256Managed();
            int numHashBytes = shaM.HashSize >> 3;
            byte[] fileHash = new byte[numHashBytes];
            _baseFileStream.Read(fileHash, 0, numHashBytes);

            // decrypt the data
            using (AesManaged aesAlg = new AesManaged())
            {
                aesAlg.Padding = PaddingMode.PKCS7;
                // convert the key to bytes
                DeriveBytes deriveBytes = new Rfc2898DeriveBytes(_key, Encoding.ASCII.GetBytes(Salt));

                aesAlg.Key = deriveBytes.GetBytes(aesAlg.KeySize >> 3);
                aesAlg.IV = deriveBytes.GetBytes(aesAlg.BlockSize >> 3);

                ICryptoTransform decryptor = aesAlg.CreateDecryptor(aesAlg.Key, aesAlg.IV);

                const int bufferSize = 1048576;
                byte[] buffer = new byte[bufferSize];

                // reset the memory stream position
                _baseMemoryStream.Position = 0;

                using (CryptoStream csDecrypt = new CryptoStream(_baseFileStream, decryptor, CryptoStreamMode.Read))
                {
                    while (true)
                    {
                        // read the next buffer of bytes
                        int numBytesRead = csDecrypt.Read(buffer, 0, bufferSize);
                        if (numBytesRead == 0) break;

                        // write those bytes to the memory stream
                        _baseMemoryStream.Write(buffer, (int) _baseMemoryStream.Position, numBytesRead);
                    }
                }

                // create a hash of the plain text and write it to the file stream
                byte[] decryptHash = shaM.ComputeHash(_baseMemoryStream.GetBuffer(), 0, (int) _baseMemoryStream.Position);

                if (!AreHashesEqual(fileHash, decryptHash))
                {
                    throw new CryptographicException(
                        string.Format(
                            "ERROR: The contents in the encrypted file ({0}) have been modified since they were created.",
                            _filePath));
                }

                // reset the memory stream position
                _baseMemoryStream.Position = 0;
            }
        }

        /// <summary>
        ///     creates an encrypted buffer from the memory stream and writes it to the file stream
        /// </summary>
        private void Encrypt()
        {
            // create a hash of the plain text and write it to the file stream
            SHA256 shaM = new SHA256Managed();
            byte[] hash = shaM.ComputeHash(_baseMemoryStream.GetBuffer(), 0, (int) _baseMemoryStream.Position);
            _baseFileStream.Write(hash, 0, hash.Length);

            // encrypt the data
            using (AesManaged aesAlg = new AesManaged())
            {
                aesAlg.Padding = PaddingMode.PKCS7;
                // convert the key to bytes
                DeriveBytes deriveBytes = new Rfc2898DeriveBytes(_key, Encoding.ASCII.GetBytes(Salt));

                aesAlg.Key = deriveBytes.GetBytes(aesAlg.KeySize >> 3);
                aesAlg.IV = deriveBytes.GetBytes(aesAlg.BlockSize >> 3);

                ICryptoTransform encryptor = aesAlg.CreateEncryptor(aesAlg.Key, aesAlg.IV);

                using (CryptoStream csEncrypt = new CryptoStream(_baseFileStream, encryptor, CryptoStreamMode.Write))
                {
                    csEncrypt.Write(_baseMemoryStream.GetBuffer(), 0, (int) _baseMemoryStream.Position);
                }
            }
        }

        /// <summary>
        ///     Clears buffers for this stream and causes any buffered data to be written to the file.
        /// </summary>
        public override void Flush()
        {
            _baseMemoryStream.Flush();
        }

        /// <summary>
        ///     opens the encrypted file stream
        /// </summary>
        private void Open(string filePath, FileAccess fileAccess)
        {
            // sanity check: we don't support ReadWrite
            if (fileAccess == FileAccess.ReadWrite)
            {
                throw new ApplicationException(
                    string.Format(
                        "ERROR: An encrypted filestream was created with an unsupported file access mode: {0}",
                        fileAccess));
            }

            // set the basic file stream parameters
            _isWriting = (fileAccess == FileAccess.Write);

            // sanity check: make sure the file exists if reading
            if (!_isWriting && !File.Exists(filePath))
            {
                throw new ApplicationException(
                    string.Format("ERROR: Could not open {0} for reading. The file does not exist.", filePath));
            }

            // open the file stream
            _baseFileStream = new FileStream(filePath, (_isWriting ? FileMode.Create : FileMode.Open));

            // read all of the contents of the file and decrypt it
            if (!_isWriting) Decrypt();
        }

        /// <summary>
        ///     Reads a block of bytes from the stream and writes the data in a given buffer.
        /// </summary>
        public override int Read(byte[] buffer, int offset, int count)
        {
            // sanity check
            if (_isWriting) throw new NotSupportedException();

            // read from the internal memory stream
            return _baseMemoryStream.Read(buffer, offset, count);
        }

        /// <summary>
        ///     Sets the current position of this stream to the given value.
        /// </summary>
        public override long Seek(long offset, SeekOrigin origin)
        {
            return _baseMemoryStream.Seek(offset, origin);
        }

        /// <summary>
        ///     Sets the length of this stream to the given value.
        /// </summary>
        public override void SetLength(long value)
        {
            _baseMemoryStream.SetLength(value);
        }

        /// <summary>
        ///     Writes a block of bytes to this stream using data from a buffer.
        /// </summary>
        public override void Write(byte[] buffer, int offset, int count)
        {
            // sanity check
            if (!_isWriting) throw new NotSupportedException();

            _baseMemoryStream.Write(buffer, offset, count);
        }

        #endregion
    }
}