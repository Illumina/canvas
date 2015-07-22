using System;
using System.IO;
using System.Xml.Serialization;
using Illumina.Common;

namespace SequencingFiles
{
    /// <summary>
    ///     Base class which provides simple serialization and deserialzation routines commonly
    ///     used for things like config files, etc.
    /// </summary>
    /// <typeparam name="TTypeToSerialize">Pass in the type of your derived class</typeparam>
    [Serializable]
    public class SerializableToFile<TTypeToSerialize> where TTypeToSerialize : SerializableToFile<TTypeToSerialize>
    {
        /// <summary>
        ///     Used when you want to deserialize this class out of a previously initialized XML-serialized file.
        ///     Derived classes are expected to provide their own public static wrapper for
        ///     this routine which simply hides the description field
        /// </summary>
        /// <param name="path">File to read and serialize into an object.</param>
        /// <param name="description">String to use in error messages referring to the object which failed serialization.</param>
        /// <returns>Deserialized object loaded from file "path"</returns>
        public static TTypeToSerialize Load(string path, string description)
        {
            return Load(path, description, null);
        }

        public static TTypeToSerialize Load(string path, string description, XmlAttributeOverrides overrides)
        {
            try
            {
                using (FileStream fs = new FileStream(path, FileMode.Open, FileAccess.Read))
                {
                    return Load(fs, overrides);
                }
            }
            catch (Exception ex)
            {
                throw AdditionalContextException.Wrap(string.Format(
                    "While loading {0} \"{1}\"", description, path),
                                                     ex
                    );
            }
        }

        /// <summary>
        ///     Used when you want to deserialize this class out of a previously initialized memory buffer.
        /// </summary>
        /// <param name="inputStream">Buffer that already contains the serialized object you want to deserialize into an object</param>
        /// <param name="overrides">A System.Xml.Serialization.XmlAttributeOverrides</param>
        /// <returns>The fully deserialized object</returns>
        public static TTypeToSerialize Load(Stream inputStream, XmlAttributeOverrides overrides)
        {
            XmlSerializer serializer = new XmlSerializer(typeof (TTypeToSerialize), overrides);
            return serializer.Deserialize(inputStream) as TTypeToSerialize;
        }


        /// <summary>
        ///     Used when you want to serialize this class into a file as XML.
        ///     Derived classes are expected to provide their own public static wrapper for
        ///     this routine which simply hides the description field.
        /// </summary>
        /// <param name="filepath">File to write to when serializing into an object.</param>
        /// <param name="description">String to use in error messages referring to the object which failed serialization.</param>
        public void Save(string filepath, string description)
        {
            try
            {
                FileInfo fi = new FileInfo(filepath);
                if (!fi.Directory.Exists)
                    fi.Directory.Create();

                if (File.Exists(filepath))
                {
                    // Eliminate possiblity of access error on fresh installs - file may come with read-only 
                    FileAttributes attributes = File.GetAttributes(filepath);
                    attributes &= ~FileAttributes.ReadOnly;
                    File.SetAttributes(filepath, attributes);
                }

                using (FileStream sw = new FileStream(filepath, FileMode.Create, FileAccess.Write))
                {
                    Save(sw);
                }
            }
            catch (Exception ex)
            {
                throw AdditionalContextException.Wrap(string.Format(
                    "While saving {0} \"{1}\"", description, filepath),
                                                     ex
                    );
            }
        }

        /// <summary>
        ///     Used when you want to serialize this class into a memory buffer.
        /// </summary>
        /// <param name="outputStream">Memory buffer to receive the serialized object</param>
        public void Save(Stream outputStream)
        {
            XmlSerializer serializer = new XmlSerializer(typeof (TTypeToSerialize));
            serializer.Serialize(outputStream, this);
        }
    }


    /// <summary>
    ///     Base class which extends SerializableToFile&lt;&gt; by adding the "Name" property,
    ///     and enforcing that this property is initialized from the filename when loaded and saved.
    /// </summary>
    /// <typeparam name="TTypeToSerialize">
    ///     Pass in the type of your derived class.  NOTE: This type must be derived from SerializableToFileWithName. Although seemingly redundant, this generic argument is still needed
    ///     in order to get a static Load() function to correctly serialize a class because a static function
    ///     has no "this" object for which to call GetType() on.
    /// </typeparam>
    [Serializable]
    public class SerializableToFileWithName<TTypeToSerialize> : SerializableToFile<TTypeToSerialize>
        where TTypeToSerialize : SerializableToFileWithName<TTypeToSerialize>
    {
        #region Events

        /// <summary>
        ///     Event that fires whenever a scan setting is added/removed or modified.
        /// </summary>
        public event EventHandler<ChangedEventArgs> Changed;

        /// <summary>
        ///     Function which fires the Changed event.
        /// </summary>
        protected void OnChanged(object sender, ChangedEventArgs args)
        {
            if (Changed != null)
                Changed(sender, args);
        }

        /// <summary>
        ///     The event arguments passed to the Changed event
        /// </summary>
        public class ChangedEventArgs : EventArgs
        {
        }

        #endregion

        private string _name;

        /// <summary>
        ///     Returns the cosmetic name for this object.  "Name" is always derived from the filename, minus extension,
        ///     that this object was originally serialized from or to.
        /// </summary>
        [XmlIgnore]
        public string Name
        {
            get { return _name; }
            set
            {
                _name = value;
                OnChanged(this, new ChangedEventArgs());
            }
        }


        /// <summary>
        ///     Used when you want to deserialize this class out of a previously initialized XML-serialized file.
        ///     Derived classes are expected to provide their own public static wrapper for
        ///     this routine which simply hides the description field
        /// </summary>
        /// <param name="path">File to read and serialize into an object.</param>
        /// <param name="description">String to use in error messages referring to the object which failed serialization.</param>
        /// <returns>Deserialized object loaded from file "path"</returns>
        public new static TTypeToSerialize Load(string path, string description)
        {
            SerializableToFileWithName<TTypeToSerialize> obj = SerializableToFile<TTypeToSerialize>.Load(path, description);

            // This assumes that typeToSerialize is really derived from this class.
            obj.Name = Path.GetFileNameWithoutExtension(path);
            return (TTypeToSerialize) obj;
        }

        /// <summary>
        ///     Used when you want to serialize this class into a file as XML.
        ///     Derived classes are expected to provide their own public static wrapper for
        ///     this routine which simply hides the description field.
        /// </summary>
        /// <param name="path">File to write to when serializing into an object.</param>
        /// <param name="description">String to use in error messages referring to the object which failed serialization.</param>
        public new void Save(string path, string description)
        {
            Name = Path.GetFileNameWithoutExtension(path);
            base.Save(path, description);
        }
    }
}