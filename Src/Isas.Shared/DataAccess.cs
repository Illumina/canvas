using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace Isas.Shared
{
    public enum DataAccessFileType
    {
        Custom,
        Position,
        Sample
    };

    public enum DataAccessPositionColumns
    {
        Position = 0,
        Cover,
        Error,
        Qscore
    };

    [DataContract]
    public class MetaDataItem
    {
        public MetaDataItem(string n, string v)
        {
            Name = n;
            Value = v;
        }

        [DataMember]
        public string Name { set; get; }

        [DataMember]
        public string Value { set; get; }
    }

    [DataContract]
    public class LogItem
    {
        public LogItem(string d, string T, string m)
        {
            Date = d;
            Time = T;
            Message = m;
        }

        [DataMember]
        public string Date { set; get; }

        [DataMember]
        public string Time { set; get; }

        [DataMember]
        public string Message { set; get; }
    }

    /// <summary>
    ///     Defines a float binary file format for saving/retrieving subsets of data
    ///     All data are stored as floats
    ///     Format is DAF[Version][NumCols][NumRows][Col1...][Col2...]
    ///     "X" values, if present, are assumed to be the first column and should be sorted smallest to largest
    /// </summary>
    public class DataAccessFile
    {
        private const int ColumnNameLength = 8; // truncate column names to this many bytes
        private static readonly char[] HeaderID = new[] { 'D', 'A', 'F' }; // first 3 bytes id this as a DataAccessFile
        private readonly string _fileName;
        private List<string> _columns = new List<string>();
        private int _numColumns;
        private int _numRows;

        private DataAccessFileType _type;
        private int _version = 1;

        public DataAccessFile(string fileName, DataAccessFileType fileType)
        {
            _fileName = fileName;
            _type = fileType;

            switch (_type)
            {
                case DataAccessFileType.Position:
                    _columns = Enum.GetNames(typeof(DataAccessPositionColumns)).ToList();
                    break;
                case DataAccessFileType.Sample:
                    break;
            }
        }

        /// <summary>
        ///     Constructs a data access file with user supplied columns
        /// </summary>
        public DataAccessFile(string fileName, List<string> userSuppliedColumns)
        {
            _fileName = fileName;
            _type = DataAccessFileType.Custom;

            _columns = userSuppliedColumns;
        }

        /// <summary>
        ///     Column names are either specified separately or read from file
        /// </summary>
        public DataAccessFile(string fileName)
        {
            _fileName = fileName;
            _type = DataAccessFileType.Custom;
        }

        public List<string> Columns
        {
            get { return _columns; }
            set { _columns = value; }
        }

        public void SaveData(float[][] data)
        {
            // make sure path exists
            string path = Path.GetDirectoryName(_fileName);
            if (!Directory.Exists(path))
            {
                Directory.CreateDirectory(path);
            }

            // do some checks
            if (data == null || data.Length == 0) return;
            _numColumns = data.Length;
            _numRows = data[0].Length;

            if (_numRows == 0) return;

            if (_numColumns != _columns.Count)
            {
                throw new ApplicationException("DataAccessFile: Bad number of columns");
            }

            using (FileStream myStream = new FileStream(_fileName, FileMode.Create, FileAccess.Write, FileShare.Read))
            {
                using (BinaryWriter writer = new BinaryWriter(myStream))
                {
                    writer.Write(_version);
                    writer.Write(HeaderID);
                    writer.Write(_numColumns);
                    foreach (string columnName in _columns)
                    {
                        // use up to 8 bytes to identify column
                        byte[] buffer = new byte[ColumnNameLength];
                        for (int charIndex = 0; charIndex < ColumnNameLength; charIndex++)
                        {
                            if (charIndex < columnName.Length) buffer[charIndex] = (byte)columnName[charIndex];
                        }
                        writer.Write(buffer);
                    }
                    writer.Write(_numRows);
                    foreach (float[] col in data)
                    {
                        if (col == null || col.Length != _numRows)
                        {
                            throw new ApplicationException("DataAccessFile: Malformed data");
                        }
                        byte[] tbuf = new byte[col.Length * sizeof(float)];
                        Buffer.BlockCopy(col, 0, tbuf, 0, col.Length * sizeof(float));
                        writer.Write(tbuf);
                    }
                }
            }
        }

        /// <summary>
        ///     Returns all columns and rows of a data access file
        /// </summary>
        /// <returns></returns>
        public float[][] GetData()
        {
            float[][] returnValues = null;

            List<string> foundColumns = new List<string>();

            using (FileStream myStream = new FileStream(_fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (BinaryReader reader = new BinaryReader(myStream))
                {
                    ReadHeader(reader, foundColumns);

                    returnValues = new float[_numColumns][];

                    // now read in all the data
                    for (int columnIndex = 0; columnIndex < _numColumns; columnIndex++)
                    {
                        returnValues[columnIndex] = new float[_numRows];
                        byte[] tbuf = reader.ReadBytes(_numRows * sizeof(float));
                        if (tbuf != null && tbuf.Length == _numRows * sizeof(float))
                        {
                            Buffer.BlockCopy(tbuf, 0, returnValues[columnIndex], 0, _numRows * sizeof(float));
                        }
                    }
                }
            }

            _columns = foundColumns;

            return returnValues;
        }

        private void ReadHeader(BinaryReader reader, List<string> foundColumns)
        {
            _version = reader.ReadInt32();
            byte[] headerCheck = reader.ReadBytes(HeaderID.Length);
            for (int headerIndex = 0; headerIndex < HeaderID.Length; headerIndex++)
            {
                if (headerCheck[headerIndex] != HeaderID[headerIndex])
                {
                    throw new ApplicationException("DataAccessFile: Not a DataAccessFile");
                }
            }

            _numColumns = reader.ReadInt32();
            if (_columns.Count > 0 && _numColumns != _columns.Count)
            {
                throw new ApplicationException("DataAccessFile: Bad number of columns");
            }

            for (int columnIndex = 0; columnIndex < _numColumns; columnIndex++)
            {
                byte[] buffer = reader.ReadBytes(ColumnNameLength);
                string colName = Encoding.ASCII.GetString(buffer);
                colName = colName.TrimEnd(new[] { '\0' }); // take off zero padding
                if (_columns.Count > 0 && _columns[columnIndex] != colName)
                {
                    throw new ApplicationException("Data Access: Bad Column Name");
                }
                foundColumns.Add(colName);
            }
            _numRows = reader.ReadInt32();
        }

        /// <summary>
        ///     Returns data from a list of specified columns in the file
        ///     By passing in integers it just takes those rows in the file
        ///     If the number of rows exceed MaxPoints the data is subsampled
        ///     At least 2 columns are always returned; the first column is the positions
        ///     StartX and EndX are 0-based
        /// </summary>
        public float[][] GetData(List<int> columns, int startX, int endX, int maxPoints, bool subsampleUsingPeakValues)
        {
            float[][] returnValues = null;
            List<string> foundColumns = new List<string>();
            if (startX < 0) startX = 0;
            if (!File.Exists(_fileName))
            {
                _columns = foundColumns;
                return returnValues;
            }
            try
            {
                using (FileStream myStream = new FileStream(_fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
                using (BinaryReader reader = new BinaryReader(myStream))
                {
                    ReadHeader(reader, foundColumns);

                    long headerPos = reader.BaseStream.Position;

                    // discover bin size
                    // read the Min and Max position in the Position column
                    long seekPos = ((int)(DataAccessPositionColumns.Position) * _numRows) * sizeof(float) + headerPos;
                    reader.BaseStream.Seek(seekPos, SeekOrigin.Begin);
                    float firstX = reader.ReadSingle();
                    float secondX = reader.ReadSingle();

                    int binSize = (int)(secondX - firstX);

                    // now adjust start/stop positions
                    if (binSize > 0)
                    {
                        startX /= binSize;
                        endX /= binSize;
                    }

                    if (endX > _numRows) endX = _numRows;
                    if (endX == 0) endX = _numRows; // convention
                    int subSample = (endX - startX + 1) / maxPoints;
                    int numberOfPoints = (endX - startX + 1) / (subSample + 1);
                    if (numberOfPoints <= 0) numberOfPoints = 1;
                    if (numberOfPoints > _numRows) numberOfPoints = _numRows;

                    if (columns != null && columns.Count > 0)
                    {
                        returnValues = new float[columns.Count][];
                        for (int columnIndex = 0; columnIndex < columns.Count; columnIndex++)
                        {
                            reader.BaseStream.Seek((columns[columnIndex] * _numRows + startX) * sizeof(float) + headerPos,
                                                   SeekOrigin.Begin);
                            returnValues[columnIndex] = new float[numberOfPoints];
                            if (subSample == 0) // return whole block
                            {
                                byte[] tbuf = reader.ReadBytes(numberOfPoints * sizeof(float));
                                if (tbuf != null && tbuf.Length >= numberOfPoints * sizeof(float))
                                    Buffer.BlockCopy(tbuf, 0, returnValues[columnIndex], 0,
                                                     numberOfPoints * sizeof(float));
                            }
                            else
                            {
                                // Subsample points.  Take the mean value (or peak value) for points in the interval
                                int bufPos = 0;
                                float peakVal = float.MinValue;

                                int pointCount = 0;
                                float total = 0;
                                for (int rowIndex = 0; rowIndex < endX - startX; rowIndex++)
                                {
                                    float tval = reader.ReadSingle();
                                    pointCount++;
                                    total += tval;
                                    peakVal = Math.Max(tval, peakVal);

                                    if (rowIndex % (subSample + 1) == 0 && bufPos < numberOfPoints)
                                    {
                                        if (subsampleUsingPeakValues)
                                        {
                                            returnValues[columnIndex][bufPos] = peakVal;
                                            peakVal = float.MinValue;
                                        }
                                        else
                                        {
                                            returnValues[columnIndex][bufPos] = total / pointCount;
                                            total = 0;
                                            pointCount = 0;
                                        }
                                        bufPos++;
                                    }
                                }
                            }
                        } // column loop
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("DataAccess Error: " + ex.Message + " :" + ex.StackTrace);
            }

            _columns = foundColumns;
            return returnValues;
        }

        /// <summary>
        ///     Returns the actual Max X value (which is not the same as the number of rows when data is binned)
        ///     This is the 0-based max, which is ChromosomeLength - 1
        /// </summary>
        public int GetMaxPosition()
        {
            int maxPosition = 0;
            if (File.Exists(_fileName))
            {
                try
                {
                    using (
                        FileStream myStream = new FileStream(_fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
                    using (BinaryReader reader = new BinaryReader(myStream))
                    {
                        ReadHeader(reader, new List<string>());
                        long headerPos = reader.BaseStream.Position;
                        // seek to next-to-last point in Position column
                        long seekPos = (((int)(DataAccessPositionColumns.Position) + 1) * _numRows - 1) * sizeof(float) +
                                       headerPos;
                        reader.BaseStream.Seek(seekPos, SeekOrigin.Begin);
                        float lastX = reader.ReadSingle();
                        maxPosition = (int)lastX;
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine("DataAccess Error: " + ex.Message + " :" + ex.StackTrace);
                }
            }
            return maxPosition;
        }
    }
}