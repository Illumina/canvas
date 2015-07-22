using System;
using System.Runtime.Serialization;

namespace Illumina.Common
{
	/// <summary>
	/// This class is used when wrapping a previous exception and rethrowing 
	/// to add additional context information to the original error.  To use, 
	/// call the static Wrap() method.  
	/// </summary>
	[Serializable]
	public class AdditionalContextException : ApplicationException
	{
		/// <summary>
		/// Use this method to wrap and re-throw exceptions with additional information about the circumstances of the original exception.  
		/// This method uses reflection to construct an exception of the same type as the original if possible
		/// so any catch statements will still receive the expected exception type.
		/// This method will create a new outer exception of the same
		/// type as the innerException argument.  The innerException argument is retained as the InnerException
		/// of the new exception.  If a new exception of the same type cannot be created, 
		/// then a new exception of type AdditionalContextException is created
		/// </summary>
		public static Exception Wrap(string additionalContextMessage, Exception innerException)
		{
			var ctor = (innerException == null)
				? null
				: innerException.GetType().GetConstructor(new [] { typeof(string), typeof(Exception) });

			var newEx = (ctor == null)
				? null
				: (Exception)ctor.Invoke(new object[] { additionalContextMessage, innerException });

			return newEx ?? new AdditionalContextException(additionalContextMessage, innerException);
		}

		private AdditionalContextException(SerializationInfo info, StreamingContext context) : base(info, context) { }

		private AdditionalContextException(string additionalContextMessage, Exception innerException)
			: base(additionalContextMessage, innerException) { }
	}
}