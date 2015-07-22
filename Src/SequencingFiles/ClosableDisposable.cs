using System;

namespace SequencingFiles
{
	interface IClosableDisposable : IDisposable
	{
		void Close();
	}

    /// <summary>
    ///     Base for classes which need to have a Close() call
    ///     triggered at dispose.
    /// </summary>
	public abstract class ClosableDisposable : IClosableDisposable
    {
        private bool _isDisposed;

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~ClosableDisposable()
        {
            Dispose(false);
        }

        /// <summary>
        ///     Closes the file
        /// </summary>
        public abstract void Close();

        // Implement IDisposable
        protected virtual void Dispose(bool disposing)
        {
            lock (this)
            {
                if (!_isDisposed)
                {
                    _isDisposed = true;
                    Close();
                }
            }
        }

        // Implement IDisposable
    }
}