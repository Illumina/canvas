using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
	public class TileSelection
	{
		private readonly Dictionary<int, HashSet<string>> _tilesPerLane = new Dictionary<int, HashSet<string>>();
		private readonly Dictionary<int, bool> _isLaneLimited = new Dictionary<int, bool>();

		private TileSelection(Dictionary<int, HashSet<string>> tilesPerLane, Dictionary<int, bool> isLaneLimited)
		{
			_tilesPerLane = tilesPerLane;
			_isLaneLimited = isLaneLimited;
		}

		public bool IsLaneLimited(int lane)
		{
			if (!_isLaneLimited.ContainsKey(lane)) throw new ArgumentException(string.Format("Tiles have not been set for this lane {0}", lane));
			return _isLaneLimited[lane];
		}

		public List<int> Lanes { get { return _tilesPerLane.Keys.ToList(); } }

		public List<string> GetTilesForLane(int lane)
		{
			if (!_tilesPerLane.ContainsKey(lane)) throw new ArgumentException(string.Format("This RunFolder does not have lane {0}", lane));
			return _tilesPerLane[lane].ToList();
		}

		public class Builder
		{
			private readonly Dictionary<int, HashSet<string>> _tilesPerLane = new Dictionary<int, HashSet<string>>();
			private readonly Dictionary<int, bool> _isLaneLimited = new Dictionary<int, bool>();
			public List<int> Lanes { get { return _tilesPerLane.Keys.ToList(); } }

			public void AddTile(int lane, string tile)
			{
				if (!_tilesPerLane.ContainsKey(lane))
					_tilesPerLane[lane] = new HashSet<string>();
				_tilesPerLane[lane].Add(tile);
				_isLaneLimited[lane] = false;
			}

			public void ExcludeTilesForLane(int lane, List<string> tiles)
			{
				if (!_tilesPerLane.ContainsKey(lane)) throw new ArgumentException(string.Format("Tiles have not been set for this lane {0}", lane));
				if (!tiles.Any()) return;
				if (!_tilesPerLane.ContainsKey(lane)) return;
				HashSet<string> limitedTiles = _tilesPerLane[lane];
				if (!limitedTiles.Overlaps(tiles)) return;
				limitedTiles.ExceptWith(tiles);
				_tilesPerLane[lane] = limitedTiles;
				_isLaneLimited[lane] = true;
			}

			public TileSelection Create()
			{
				return new TileSelection(_tilesPerLane, _isLaneLimited);
			}
		}
	}
}
