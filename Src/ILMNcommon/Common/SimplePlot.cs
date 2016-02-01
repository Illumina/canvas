using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.IO;
using System.Xml;
using System.Xml.Linq;

// Disable "Comparison made to same variable" warning, since we use it as a much faster check for NaN:
#pragma warning disable 1718

namespace Illumina.Common
{
    public class Padding
    {
        public int BottomPadding;
        public int LeftPadding;
        public int RightPadding;
        public int TopPadding;
    }

    public class SimplePlot
    {
        #region Members

        public static object PlotLock = new object();
        public int BottomPadding = 39;
        protected int Height = 100;
        public int LeftPadding = 65;
        public int MinBottomPadding = 10;
        public int RightPadding = 40; //25;
        public int TopPadding = 30;
        protected int Width = 100;
        public float XAxisOffset = 0;
        private bool _drawRightYAxisLabels = false;
        protected int ProtectedBottomPadding = 39;
        protected int ProtectedLeftPadding = 40;
        protected float PixelPerX = 1.0f;
        protected float PixelPerY = 1.0f;
        protected int ProtectedRightPadding = 40; //25;
        protected int ProtectedTopPadding = 30;
        public bool ShowTickMarksX { get; set; }
        public bool ShowTickMarksY { get; set; }

        #endregion

        private readonly List<GraphObject> _graphObjects = new List<GraphObject>();
        private float _xAxisLabelRotation;
        public List<AxisLabel> XAxisLabels = new List<AxisLabel>();
        public List<AxisLabel> YAxisLabels = new List<AxisLabel>();
        private bool _showTickLabels = true;
        private string _yUnits = "";
        private Color _axesLabelColor = Color.Black;
        private Font _axesLabelFont = new Font("Microsoft Sans Serif", 9);
        private Color _backColor = Color.White;
        private Font _captionFont = new Font("Microsoft Sans Serif", 14);
        private Color _captionLabelColor = Color.Red;
        private Color _foreColor = Color.Black;
        private Color _gridColor = Color.LightGray;
        private bool _hideXLabels;
        private Color _labelColor = Color.Black;
        private Font _labelFont = new Font("Microsoft Sans Serif", 7);
        private float _maxX = 100;
        private float _maxY = 100;
        private bool _showGridLines = true;
        private Color _titleColor = Color.Black;
        private Font _titleFont = new Font("Microsoft Sans Serif", 9);

        public SimplePlot(int width, int height)
        {
            ShowTickMarksX = true;
            ShowTickMarksY = true;
            SetSize(width, height);
        }

        public Padding Padding
        {
            get
            {
                Padding padding = new Padding
                                      {
                                          LeftPadding = ProtectedLeftPadding,
                                          RightPadding = ProtectedRightPadding,
                                          BottomPadding = ProtectedBottomPadding,
                                          TopPadding = ProtectedTopPadding
                                      };
                return padding;
            }
        }

        public string Title { get; set; }

        public string YLabel { get; set; }

        public string XLabel { get; set; }

        public string YUnits
        {
            get { return _yUnits; }
            set { _yUnits = value; }
        }

        public bool ShowGridLines
        {
            get { return _showGridLines; }
            set { _showGridLines = value; }
        }

        public bool HideXLabels
        {
            get { return _hideXLabels; }
            set { _hideXLabels = value; }
        }

        public float MinX { get; set; }

        public float MinY { get; set; }

        public float MaxX
        {
            get { return _maxX; }
            set { _maxX = value; }
        }

        public float MaxY
        {
            get { return _maxY; }
            set { _maxY = value; }
        }

        public bool ShowXTickLabels
        {
            get { return _showTickLabels; }
            set { _showTickLabels = value; }
        }

        public Color GridColor
        {
            get { return _gridColor; }
            set { _gridColor = value; }
        }

        public Color BackColor
        {
            get { return _backColor; }
            set { _backColor = value; }
        }

        public Color ForeColor
        {
            get { return _foreColor; }
            set
            {
                _foreColor = value;
                LabelColor = TitleColor = AxesLabelColor = value;
            }
        }

        public Color LabelColor
        {
            get { return _labelColor; }
            set { _labelColor = value; }
        }

        public Color TitleColor
        {
            get { return _titleColor; }
            set { _titleColor = value; }
        }

        public Color AxesLabelColor
        {
            get { return _axesLabelColor; }
            set { _axesLabelColor = value; }
        }

        public Color CaptionLabelColor
        {
            get { return _captionLabelColor; }
            set { _captionLabelColor = value; }
        }

        public Font LabelFont
        {
            get { return _labelFont; }
            set { _labelFont = value; }
        }

        public Font TitleFont
        {
            get { return _titleFont; }
            set { _titleFont = value; }
        }

        public Font AxesLabelFont
        {
            get { return _axesLabelFont; }
            set { _axesLabelFont = value; }
        }

        public Font CaptionFont
        {
            get { return _captionFont; }
            set { _captionFont = value; }
        }

        public float XAxisLabelRotation
        {
            get { return _xAxisLabelRotation; }
            set { _xAxisLabelRotation = value; }
        }

        public float GetPixelsPerX()
        {
            return PixelPerX;
        }

        public void SetSize(int width, int height)
        {
            Width = width;
            Height = height;
            if (OKToRender())
            {
                RecalculateAxes();
                RecalculateObjects();
            }
        }

        protected static float FindMin(float absmin, float tickspacing)
        {
            // this will point to the first visible tickblock (since absmin may not be on a tick position)
            float firstVisibleTickMin = 0;
            if (tickspacing <= 0) return 0;
            // number of tickblocks before absmin
            int tickBlocks = (int) (Math.Abs(absmin)/tickspacing);

            if (absmin > 0)
            {
                firstVisibleTickMin = (tickBlocks + 1)*tickspacing;
            }
            else
            {
                firstVisibleTickMin = -tickBlocks*tickspacing;
            }
            return firstVisibleTickMin;
        }

        protected static float GetLogTickSpacingEstimate(float minv, float maxv)
        {
            float diff = (float) (Math.Log10(maxv) - Math.Log10(minv));
            int powrMin = (int) Math.Floor(Math.Log10(diff/10.0));

            float factor = (float) (Math.Pow(10, powrMin + 1));
            if (factor >= 1) return factor;
            return 1;
        }

        protected static float GetTickSpacingEstimate(float minv, float maxv)
        {
            float diff = (maxv - minv);
            int powrMin = (int) Math.Floor(Math.Log10(diff/10.0));
            int powrMax = (int) Math.Floor(Math.Log10(diff/3.0));

            if (powrMin != powrMax)
            {
                return (float) (Math.Pow(10, powrMin + 1));
            }

            return (float) (2.0*Math.Pow(10, powrMin));
        }

        protected virtual void DrawGridLines(Graphics g, float labelSpacing)
        {
            if (!OKToRender())
                return;
            int ypos = Height - ProtectedBottomPadding;

            float ylabelspacing = GetTickSpacingEstimate(MinY, MaxY);
            SizeF size = g.MeasureString("0", LabelFont);
            while (size.Height > (ConvertYToPixel(0) - ConvertYToPixel(ylabelspacing)))
            {
                ylabelspacing = ylabelspacing*2;
            }

            if (labelSpacing == 0 || ylabelspacing == 0)
                return;

            double xval = FindMin(MinX, labelSpacing);
            float yval = FindMin(MinY, ylabelspacing);

            double oldxval = xval;
            float oldyval = yval;

            using (Pen gridPen = new Pen(new SolidBrush(GridColor), 1))
            {
                using (Pen forePen = new Pen(new SolidBrush(ForeColor), 1))
                {
                    // draw the x-grid lines
                    float tMax = MaxX;
                    int xpos;
                    if (ShowGridLines)
                    {
                        while (xval < tMax)
                        {
                            xpos = (int) ConvertXToPixel(xval);

                            g.DrawLine(Math.Abs(xval) < labelSpacing/2.0f ? forePen : gridPen, xpos, ProtectedTopPadding, xpos,
                                       ypos);
                            xval += labelSpacing;
                        }
                    }

                    if (ShowTickMarksX)
                    {
                        xval = oldxval;
                        while (xval < tMax)
                        {
                            xpos = (int) ConvertXToPixel(xval);
                            g.DrawLine(forePen, xpos, ypos - 4, xpos, ypos);
                            xval += labelSpacing;
                        }
                    }

                    if (ShowGridLines)
                    {
                        // draw the y-grid lines
                        while (yval < MaxY)
                        {
                            ypos = (int) ConvertYToPixel(yval);
                            g.DrawLine(Math.Abs(yval) < ylabelspacing/2.0f ? forePen : gridPen, ProtectedLeftPadding, ypos,
                                       Width - ProtectedRightPadding, ypos);
                            yval += ylabelspacing;
                        }
                    }

                    if (ShowTickMarksY)
                    {
                        yval = oldyval;
                        while (yval < MaxY)
                        {
                            ypos = (int) ConvertYToPixel(yval);
                            g.DrawLine(forePen, ProtectedLeftPadding, ypos, ProtectedLeftPadding + 4, ypos);
                            g.DrawLine(forePen, Width - ProtectedRightPadding - 4, ypos, Width - ProtectedRightPadding, ypos);
                            yval += ylabelspacing;
                        }
                    }
                }
            }
        }

        protected virtual void PlotLabels(Graphics g, float labelspacing)
        {
            if (!OKToRender())
                return;
            int ypos = Height - ProtectedBottomPadding; //y2Pxl(0.0);

            float ylabelspacing = GetTickSpacingEstimate(MinY, MaxY);

            SizeF size = g.MeasureString("0", LabelFont);
            while (size.Height > (ConvertYToPixel(0) - ConvertYToPixel(ylabelspacing)))
            {
                ylabelspacing = ylabelspacing*2;
            }

            if (labelspacing == 0 || ylabelspacing == 0)
                return;
            double xval = FindMin(MinX, labelspacing);
            float yval = FindMin(MinY, ylabelspacing);

            double oldxval = xval;
            float oldyval = yval;

            // draw the x grid labels
            xval = oldxval;

            StringFormat format = new StringFormat
                                      {
                                          Alignment = StringAlignment.Center,
                                          LineAlignment = StringAlignment.Center
                                      };
            using (Brush labelBrush = new SolidBrush(LabelColor))
            {
                using (Brush titleBrush = new SolidBrush(TitleColor))
                {
                    using (Brush axesLabelBrush = new SolidBrush(AxesLabelColor))
                    {
                        int xpos;
                        string cstr;
                        SizeF strsize;
                        if (_showTickLabels)
                        {
                            float tMax = MaxX;
                            while (xval <= tMax)
                            {
                                double displayXVal = xval + XAxisOffset;
                                xpos = (int) ConvertXToPixel(xval);
                                if (displayXVal == (int) displayXVal)
                                    cstr = string.Format("{0}", (int) displayXVal);
                                else
                                {
                                    if (labelspacing < 0.06 || xval < .01)
                                        cstr = string.Format("{0:F3}", displayXVal);
                                    else if (labelspacing < 0.6 || xval < .1)
                                        cstr = string.Format("{0:F2}", displayXVal);
                                    else
                                        cstr = string.Format("{0:F1}", displayXVal);
                                }
                                if (cstr.EndsWith("000000000")) cstr = cstr.Substring(0, cstr.Length - 9) + "G";
                                if (cstr.EndsWith("000000")) cstr = cstr.Substring(0, cstr.Length - 6) + "M";
                                if (cstr.EndsWith("000")) cstr = cstr.Substring(0, cstr.Length - 3) + "K";

                                strsize = g.MeasureString(cstr, LabelFont);
                                if (!_hideXLabels)
                                {
                                    g.DrawString(cstr, LabelFont, labelBrush,
                                                 new Rectangle(xpos - (int) strsize.Width/2 - 2,
                                                               Height - ProtectedBottomPadding + 1,
                                                               (int) strsize.Width + 4, (int) strsize.Height), format);
                                }
                                xval = xval + labelspacing;
                            }
                        }
                        else if (XAxisLabels.Count > 0 && _hideXLabels == false)
                        {
                            int minimumXPos = int.MinValue;
                            float sin = (float) Math.Sin(Math.PI*Math.Abs(_xAxisLabelRotation)/180);
                            foreach (AxisLabel label in XAxisLabels)
                            {
                                xpos = (int) ConvertXToPixel(label.Position);
                                int ytry = Height - ProtectedBottomPadding + 1;
                                strsize = g.MeasureString(label.Value, LabelFont);
                                Matrix saveCoordFrame1 = g.Transform.Clone();
                                int xtry;
                                int xkeepout;
                                if (_xAxisLabelRotation != 0)
                                {
                                    xtry = xpos - 2 - (int) strsize.Width;
                                    ytry = Height - ProtectedBottomPadding - 1;
                                    xkeepout = (int) (strsize.Height/sin) + 1;

                                    Matrix m1 = new Matrix();
                                    m1.RotateAt(_xAxisLabelRotation, new PointF(xpos, Height - ProtectedBottomPadding));
                                    g.MultiplyTransform(m1);

                                    if (xpos > minimumXPos)
                                    {
                                        g.DrawString(label.Value, LabelFont, labelBrush,
                                                     new Rectangle(xtry, ytry, (int) (strsize.Width + 1),
                                                                   (int) (strsize.Height + 1)), format);
                                        minimumXPos = xpos + xkeepout;
                                    }
                                }
                                else
                                {
                                    xkeepout = (int) (strsize.Width/2 + 2);
                                    xtry = xpos - (int) strsize.Width/2 - 2;
                                    if (xtry > minimumXPos)
                                    {
                                        g.DrawString(label.Value, LabelFont, labelBrush,
                                                     new Rectangle(xtry, Height - ProtectedBottomPadding + 1,
                                                                   (int) strsize.Width + 4, (int) strsize.Height),
                                                     format);
                                        minimumXPos = xpos + xkeepout;
                                    }
                                }
                                g.Transform = saveCoordFrame1;
                            }
                        }

                        // the x axis label;
                        if (!_hideXLabels)
                        {
                            g.DrawString(XLabel, AxesLabelFont, axesLabelBrush,
                                         new Rectangle(ProtectedLeftPadding, Height - ProtectedBottomPadding + 2,
                                                       Width - ProtectedLeftPadding - ProtectedRightPadding, ProtectedBottomPadding), format);
                        }

                        // draw the title
                        g.DrawString(Title, TitleFont, titleBrush,
                                     new Rectangle(ProtectedLeftPadding, 0, Width - ProtectedLeftPadding - ProtectedRightPadding, ProtectedTopPadding),
                                     format);

                        // draw the y grid labels
                        format.Alignment = StringAlignment.Far;
                        float maxwidth = 0;
                        ypos = (int) ConvertYToPixel(MaxY) - 10;
                        cstr = _yUnits;
                        strsize = g.MeasureString(cstr, LabelFont);
                        format.Alignment = StringAlignment.Far;
                        g.DrawString(cstr, LabelFont, labelBrush,
                                     new Rectangle(0, ypos - (int) strsize.Height/2 - 1, ProtectedLeftPadding - 2,
                                                   (int) strsize.Height + 2), format);
                        format.Alignment = StringAlignment.Near;
                        g.DrawString(cstr, LabelFont, labelBrush,
                                     new Rectangle(Width - ProtectedRightPadding + 2, ypos - (int) strsize.Height/2 - 1,
                                                   ProtectedRightPadding - 2, (int) strsize.Height + 2), format);
                        yval = oldyval;
                        while (yval < MaxY)
                        {
                            ypos = (int) ConvertYToPixel(yval);
                            if (yval == (int) yval)
                                cstr = string.Format("{0}", (int) yval);
                            else
                            {
                                if (ylabelspacing < 0.06)
                                    cstr = string.Format("{0:F3}", yval);
                                else if (ylabelspacing < .6)
                                    cstr = string.Format("{0:F2}", yval);
                                else
                                    cstr = string.Format("{0:F1}", yval);
                            }
                            if (cstr.EndsWith("000000")) cstr = cstr.Substring(0, cstr.Length - 6) + "M";
                            if (cstr.EndsWith("000")) cstr = cstr.Substring(0, cstr.Length - 3) + "K";
                            strsize = g.MeasureString(cstr, LabelFont);
                            if (strsize.Width > maxwidth) maxwidth = strsize.Width;
                            //Y scale
                            format.Alignment = StringAlignment.Far;
                            int top = ypos - (int) strsize.Height/2 - 1;
                            if (top < 0) top = 0;
                            g.DrawString(cstr, LabelFont, labelBrush,
                                         new Rectangle(0, top, ProtectedLeftPadding - 2, (int) strsize.Height + 2), format);
                            format.Alignment = StringAlignment.Near;
                            if (_drawRightYAxisLabels)
                                g.DrawString(cstr, LabelFont, labelBrush,
                                             new Rectangle(Width - ProtectedRightPadding + 2, top, ProtectedRightPadding - 2,
                                                           (int) strsize.Height + 2), format);
                            yval += ylabelspacing;
                        }

                        // the y axis label;
                        format.Alignment = StringAlignment.Center;
                        format.LineAlignment = StringAlignment.Center;
                        Matrix saveCoordFrame = g.Transform.Clone();
                        Matrix m = new Matrix();
                        m.Rotate(-90F);
                        g.MultiplyTransform(m);
                        g.TranslateTransform(
                            (ProtectedLeftPadding - 20 - maxwidth)*saveCoordFrame.Elements[0],
                            (Height - ProtectedBottomPadding)*saveCoordFrame.Elements[3],
                            MatrixOrder.Append);
                        g.DrawString(YLabel, AxesLabelFont, axesLabelBrush,
                                     new Rectangle(0, 0, Height - ProtectedBottomPadding - ProtectedTopPadding, 20), format);
                        g.Transform = saveCoordFrame;
                    }
                }
            }
        }

        protected virtual void DrawBackground(Graphics g)
        {
            using (Brush backBrush = new SolidBrush(BackColor))
            {
                using (Pen forePen = new Pen(ForeColor))
                {
                    g.FillRectangle(backBrush, 0, 0, Width, Height);

                    g.FillRectangle(backBrush, ProtectedLeftPadding, ProtectedTopPadding,
                                    Width - ProtectedRightPadding - ProtectedLeftPadding + 1,
                                    Height - ProtectedBottomPadding - ProtectedTopPadding + 1);

                    Region oldclip = g.Clip;
                    g.SetClip(new Rectangle(ProtectedLeftPadding, ProtectedTopPadding,
                                            Width - ProtectedRightPadding - ProtectedLeftPadding + 1,
                                            Height - ProtectedBottomPadding - ProtectedTopPadding + 1));

                    g.SetClip(oldclip, CombineMode.Replace);

                    float labelSpacing = GetRefinedTickLabelSpacing(g, LabelFont);

                    PlotLabels(g, labelSpacing);
                    DrawGridLines(g, labelSpacing);

                    g.DrawRectangle(forePen, ProtectedLeftPadding, ProtectedTopPadding,
                                    Width - ProtectedRightPadding - ProtectedLeftPadding + 1,
                                    Height - ProtectedBottomPadding - ProtectedTopPadding + 1);

                    g.SetClip(new Rectangle(ProtectedLeftPadding, ProtectedTopPadding,
                                            Width - ProtectedRightPadding - ProtectedLeftPadding + 1,
                                            Height - ProtectedBottomPadding - ProtectedTopPadding + 1));
                    foreach (GraphObject obj in _graphObjects)
                        obj.Draw(g);
                    g.SetClip(oldclip, CombineMode.Replace);
                }
            }
        }

        /// <summary>
        ///     Smart estimator of tick spacing, factoring in the font size to prevent collisions!
        /// </summary>
        private float GetRefinedTickLabelSpacing(Graphics g, Font labelFont)
        {
            float labelSpacing = GetTickSpacingEstimate(MinX, MaxX);
            while (true)
            {
                bool collision = false;
                double xval = FindMin(MinX, labelSpacing);
                int oldXPos = -99;
                while (xval <= MaxX)
                {
                    int xpos = (int) ConvertXToPixel(xval);
                    string cstr;
                    if (xval == (int) xval)
                        cstr = string.Format("{0}", (int) xval);
                    else
                    {
                        if (labelSpacing < 0.06 || xval < .01)
                            cstr = string.Format("{0:F3}", xval);
                        else if (labelSpacing < 0.6 || xval < .1)
                            cstr = string.Format("{0:F2}", xval);
                        else
                            cstr = string.Format("{0:F1}", xval);
                    }
                    if (cstr.EndsWith("000000000")) cstr = cstr.Substring(0, cstr.Length - 9) + "G";
                    if (cstr.EndsWith("000000")) cstr = cstr.Substring(0, cstr.Length - 6) + "M";
                    if (cstr.EndsWith("000")) cstr = cstr.Substring(0, cstr.Length - 3) + "K";

                    SizeF strsize = g.MeasureString(cstr, labelFont);
                    if (strsize.Width >= (xpos - oldXPos)) collision = true;
                    oldXPos = xpos;
                    xval += labelSpacing;
                }
                if (!collision) break;
                labelSpacing *= 2;
            }
            return labelSpacing;
        }

        protected virtual void SetUpXAxis(float minx, float maxx)
        {
            MinX = minx;
            MaxX = maxx;
            if (OKToRender())
            {
                RecalculateAxes();
                RecalculateObjects();
            }
        }

        protected virtual void SetUpYAxis(float miny, float maxy)
        {
            MinY = miny;
            MaxY = maxy;
            if (OKToRender())
            {
                RecalculateAxes();
                RecalculateObjects();
            }
        }

        /// <summary>
        ///     whenever the control is resized, or the axes limits are changed, this function is called
        /// </summary>
        protected virtual void RecalculateAxes()
        {
            ProtectedTopPadding = TopPadding;
            ProtectedRightPadding = RightPadding;
            ProtectedLeftPadding = LeftPadding;

            if (MaxX == MinX)
            {
                MaxX += 1;
                MinX -= 1;
            }
            if (MaxY == MinY)
            {
                MaxY += 1;
                MinY -= 1;
            }

            PixelPerX = (Width - ProtectedLeftPadding - ProtectedRightPadding)/(MaxX - MinX);
            PixelPerY = (Math.Max(1, Height - ProtectedBottomPadding - ProtectedTopPadding))/(MaxY - MinY);
        }

        private int SetBottomPadding()
        {
            int padding = BottomPadding;
            if (_hideXLabels)
            {
                return (MinBottomPadding);
            }
            if (XAxisLabels != null && _xAxisLabelRotation != 0)
            {
                using (Bitmap b = new Bitmap(10, 10))
                using (Graphics g = Graphics.FromImage(b))
                {
                    float maxLabelWidth = 0;
                    float maxLabelHeight = 0;
                    foreach (AxisLabel label in XAxisLabels)
                    {
                        SizeF strsize = g.MeasureString(label.Value, LabelFont);
                        if (strsize.Width > maxLabelWidth) maxLabelWidth = strsize.Width;
                        if (strsize.Height > maxLabelHeight) maxLabelHeight = strsize.Height;
                    }
                    padding = Math.Max(BottomPadding,
                                       (int)
                                       (2 + maxLabelHeight +
                                        maxLabelWidth*Math.Sin((Math.PI*Math.Abs(_xAxisLabelRotation)/180))));
                }
            }
            return padding;
        }

        protected virtual void RecalculateObjects()
        {
            foreach (GraphObject obj in _graphObjects)
                obj.Calculate();
        }

        public virtual void AutoScaleAxes()
        {
            float minx = float.MaxValue;
            float maxx = float.MinValue;
            float miny = float.MaxValue;
            float maxy = float.MinValue;

            foreach (GraphObject obj in _graphObjects)
            {
                float val = obj.GetMinXVal();
                if (val == float.MaxValue)
                    continue;
                if (val < minx)
                    minx = val;
                val = obj.GetMaxXVal();
                if (val > maxx)
                    maxx = val;
                val = obj.GetMinYVal();
                if (val < miny)
                    miny = val;
                val = obj.GetMaxYVal();
                if (val > maxy)
                    maxy = val;
            }

            float diff;
            if (minx == maxx)
            {
                maxx += 1;
                minx -= 1;
            }
            if (minx < maxx)
            {
                diff = maxx - minx;
                MinX = minx - 0.1f*diff;
                MaxX = maxx + 0.1f*diff;
            }
            if (miny == maxy)
            {
                maxy += 1;
                miny -= 1;
            }
            if (miny < maxy)
            {
                diff = maxy - miny;
                MinY = miny - 0.1f*diff;
                MaxY = maxy + 0.1f*diff;
            }

            if (OKToRender())
            {
                RecalculateAxes();
                RecalculateObjects();
            }
        }

        private bool OKToRender()
        {
            ProtectedBottomPadding = SetBottomPadding();
            return true;
            //if (Width > 10 && Width > leftPadding_ + rightPadding_ &&
            //    Height > 10 && Height > topPadding_ + bottomPadding_) return true;
            //else return false;
        }

        public Bitmap DrawToBitmap()
        {
            Bitmap bitmap = new Bitmap(Width, Height, PixelFormat.Format32bppArgb);
            if (OKToRender())
            {
                using (Graphics g = Graphics.FromImage(bitmap))
                    DrawBackground(g);
            }
            return bitmap;
        }

        public void DrawToBitmap(Bitmap bitmap, int x, int y)
        {
            if (OKToRender())
            {
                using (Graphics g = Graphics.FromImage(bitmap))
                {
                    g.TranslateTransform(x, y);
                    DrawBackground(g);
                    g.ResetTransform();
                }
            }
        }

        public void SetUpAxes(float minx, float maxx, float miny, float maxy)
        {
            MinX = minx;
            MaxX = maxx;
            MinY = miny;
            MaxY = maxy;
            SetUpAxes();
        }

        public void SetUpAxes()
        {
            if (OKToRender())
            {
                RecalculateAxes();
                RecalculateObjects();
            }
        }

        public float ConvertXToPixel(float xval)
        {
            if (xval != xval) return 0;
            return ProtectedLeftPadding + (int) ((xval - MinX)*PixelPerX);
        }

        public float ConvertXToPixel(int xval)
        {
            if (xval != xval) return 0;
            return ProtectedLeftPadding + (int) ((xval - (double) MinX)*PixelPerX);
        }

        public float ConvertXToPixel(double xval)
        {
            if (xval != xval) return 0;
            return ProtectedLeftPadding + (int) ((xval - MinX)*PixelPerX);
        }

        public float ConvertYToPixel(float yval)
        {
            if (yval != yval) return 0;
            return (Height - ProtectedBottomPadding) - (int) ((yval - MinY)*PixelPerY);
        }

        public float ConvertYToPixel(double yval)
        {
            if (yval != yval) return 0;
            return (Height - ProtectedBottomPadding) - (int) ((yval - MinY)*PixelPerY);
        }

        public float ConvertPixelToX(int val)
        {
            return (val - ProtectedLeftPadding)/PixelPerX + MinX;
        }

        public float ConvertPixelToY(int val)
        {
            return MaxY - (val - ProtectedBottomPadding)/PixelPerY;
        }

        public PointF GraphToScreen(PointF fpoint)
        {
            return new PointF(ConvertXToPixel(fpoint.X), ConvertYToPixel(fpoint.Y));
        }

        public PointF GraphToScreen(PointD fpoint)
        {
            return new PointF(ConvertXToPixel(fpoint.X), ConvertYToPixel(fpoint.Y));
        }

        internal void RemoveObject(GraphObject graphObject)
        {
            _graphObjects.Remove(graphObject);
        }

        internal void AddObject(GraphObject graphObject)
        {
            _graphObjects.Add(graphObject);
        }

        public void SaveToXML(string fileName, string colName)
        {
            XDocument xmlDocument = new XDocument(
                new XComment("Illumina RTA Data"));
            XElement root = new XElement("Data");
            int key = 0;
            foreach (GraphObject obj in _graphObjects)
                if (obj is SimpleBoxPlotObject)
                {
                    if (((SimpleBoxPlotObject) obj).SerializeToXML)
                    {
                        XElement row = obj.SaveToXElement(colName, (++key).ToString());
                        if (row != null) root.Add(row);
                    }
                }

            string dir = Path.GetDirectoryName(fileName);
            if (!Directory.Exists(dir)) Directory.CreateDirectory(dir);
            xmlDocument.Add(root);
            xmlDocument.Declaration = new XDeclaration("1.0", "utf-8", "true");
            using (FileStream fs = File.Open(fileName, FileMode.Create, FileAccess.Write, FileShare.None))
            {
                using (XmlWriter xw = XmlWriter.Create(fs))
                {
                    xmlDocument.Save(xw);
                    xw.Flush();
                    xw.Close();
                }
                fs.Flush();
                fs.Close();
            }
        }
    }


    public class SimpleBoxPlotObject : GraphObject
    {
        #region Members

        private readonly int _n;
        private readonly GraphPointSet _gps;
        private readonly float _iqr;
        private readonly float _max;
        private readonly float _min;
        private readonly float _p25;
        private readonly float _p50;
        private readonly float _p75;
        private readonly float _x;
        public string Label;
        private Color _color = Color.Blue;
        private Color _mediancolor = Color.Red;
        private bool _serializeToXML = true;

        #endregion

        public SimpleBoxPlotObject(SimplePlot graph, int x,
                                   List<float> vals, string label)
        {
            Label = label;
            SimplePlot = graph;
            SimplePlot.AddObject(this);
            _x = x;
            float[] valarray = vals.ToArray();
            Array.Sort(valarray);
            int n = valarray.Length;
            _n = n;
            if (n > 0)
            {
                _p25 = MathSupportFunctions.PercentilePreSorted(valarray, 25); // valarray[N / 4];
                _p50 = MathSupportFunctions.PercentilePreSorted(valarray, 50);
                _p75 = MathSupportFunctions.PercentilePreSorted(valarray, 75);
                _iqr = _p75 - _p25;
            }

            _gps = new GraphPointSet(graph);

            float llimit = _p25 - 1.5f*_iqr;
            float ulimit = _p75 + 1.5f*_iqr;

            _min = float.MaxValue;
            _max = float.MinValue;

            for (int i = 0; i < n; i++)
            {
                float v = vals[i];
                if (v < llimit || v > ulimit)
                {
                    _gps.AddPoint(new PointF(x, v), i);
                }
                else
                {
                    if (v < _min)
                        _min = v;
                    if (v > _max)
                        _max = v;
                }
            }
        }

        public Color Color
        {
            get { return _color; }
            set { _color = value; }
        }

        public Color MedianColor
        {
            get { return _mediancolor; }
            set { _mediancolor = value; }
        }

        public bool SerializeToXML
        {
            get { return _serializeToXML; }
            set { _serializeToXML = value; }
        }

        public override void Calculate()
        {
        }

        public override void Draw(Graphics g)
        {
            if (_max == float.MinValue || _min == float.MinValue) return;
            if (_max != _max || _min != _min)
                return;
            float left = SimplePlot.ConvertXToPixel(_x - 0.4f);
            float right = SimplePlot.ConvertXToPixel(_x + 0.4f);
            float center = SimplePlot.ConvertXToPixel(_x);
            float p25 = SimplePlot.ConvertYToPixel(_p25);
            float p50 = SimplePlot.ConvertYToPixel(_p50);
            float p75 = SimplePlot.ConvertYToPixel(_p75);
            float upper = SimplePlot.ConvertYToPixel(_max);
            float lower = SimplePlot.ConvertYToPixel(_min);

            SimplePlot.BottomPadding = 60;

            try
            {
                using (Pen pen = new Pen(_mediancolor))
                    g.DrawLine(pen, left, p50, right, p50);
                using (Pen pen = new Pen(_color))
                {
                    g.DrawLine(pen, center, p75, center, upper);
                    g.DrawLine(pen, (center + left)/2, upper, (center + right)/2, upper);
                    g.DrawLine(pen, center, p25, center, lower);
                    g.DrawLine(pen, (center + left)/2, lower, (center + right)/2, lower);
                    g.DrawRectangle(pen, left, p75,
                                    right - left + 1, p25 - p75 + 1);
                }
            }
            catch (OverflowException)
            {
                return;
            }
            float x = left;
            float w = right - left + 1;
            const float h = 12;
            float sw = SimplePlot.ConvertXToPixel(_x + 0.5f) - SimplePlot.ConvertXToPixel(_x - 0.5f) + 1;
            float y = SimplePlot.ConvertYToPixel(SimplePlot.MinY) - h;
            float size = 12;

            StringFormat format = new StringFormat
                                      {
                                          Alignment = StringAlignment.Center,
                                          LineAlignment = StringAlignment.Center
                                      };
            do
            {
                using (Font font = new Font(FontFamily.GenericMonospace, size))
                {
                    SizeF ssize = g.MeasureString(_n.ToString(), font);
                    if (ssize.Width <= w)
                    {
                        g.DrawString(_n.ToString(), font, Brushes.Blue, new RectangleF(x, y, w, h), format);
                        break;
                    }
                    size -= 0.5f;
                }
            } while (size > 1);

            y = SimplePlot.ConvertYToPixel(SimplePlot.MinY) + 2;
            size = 12;

            Region oldclip = g.Clip;
            g.ResetClip();
            string label2 = Label;
            if (label2.Length < 4)
                label2 = "000.";
            using (Brush brush = new SolidBrush(SimplePlot.LabelColor))
                do
                {
                    using (Font font = new Font(FontFamily.GenericMonospace, size))
                    {
                        SizeF ssize = g.MeasureString(label2, font);
                        if (ssize.Width <= sw)
                        {
                            g.DrawString(Label, font, brush, new RectangleF(x, y, sw, h), format);
                            break;
                        }
                        size -= 0.5f;
                    }
                } while (size > 1);
            g.SetClip(oldclip, CombineMode.Replace);
        }

        public override float GetMaxYVal()
        {
            float max = _gps.GetMaxYVal();
            if (_max > max)
                max = _max;
            return max;
        }

        public override float GetMaxXVal()
        {
            return _x + 1;
        }

        public override float GetMinXVal()
        {
            return _x - 1;
        }

        public override float GetMinYVal()
        {
            float min = _gps.GetMinYVal();
            if (_min < min)
                min = _min;
            return min;
        }

        public override XElement SaveToXElement(string colName, string rowKey)
        {
            if (_n > 0)
                return new XElement(colName, new XAttribute("key", rowKey), //new XAttribute("n", this._N),
                                    new XAttribute("min", _min.ToString("0.000")),
                                    new XAttribute("max", _max.ToString("0.000")),
                                    //new XAttribute("iqr", this._iqr.ToString("0.000")), new XAttribute("p25", this._p25.ToString("0.000")), new XAttribute("p75", this._p75.ToString("0.000")),
                                    new XAttribute("p50", _p50.ToString("0.000")));
            return null;
        }
    }

    public class GraphPointSet : GraphObject
    {
        protected float ProtectedBarWidth = 0;
        public Color Color = Color.Blue;
        protected float ProtectedDotWidth = 3;
        protected bool ProtectedDrawBars = false;
        protected bool ProtectedDrawDataLabels = false;
        protected ArrayList Indexes = new ArrayList();
        protected bool ProtectedLinkDots = false;
        private float _penWidth = 1;
        protected List<PointD> ProtectedPoints = new List<PointD>();
        protected List<PointF> ScreenPoints = new List<PointF>();
        public bool VerticalGradient = false;
        public Color VerticalGradientFrom = Color.FromArgb(Convert.ToInt32("FF68C930", 16));
        public Color VerticalGradientTo = Color.FromArgb(Convert.ToInt32("FF2C8E49", 16));

        public GraphPointSet()
        {
        }

        public GraphPointSet(SimplePlot simplePlot)
            : base(simplePlot)
        {
            SimplePlot.AddObject(this);
        }

        public GraphPointSet(SimplePlot simplePlot, Color color)
            : base(simplePlot)
        {
            SimplePlot = simplePlot;
            SimplePlot.AddObject(this);
            Color = color;
        }

        public List<PointD> Points
        {
            get { return ProtectedPoints; }
            set { ProtectedPoints = value; }
        }

        public float DotWidth
        {
            get { return ProtectedDotWidth; }
            set { ProtectedDotWidth = value; }
        }

        public bool LinkDots
        {
            get { return ProtectedLinkDots; }
            set { ProtectedLinkDots = value; }
        }

        public bool DrawBars
        {
            get { return ProtectedDrawBars; }
            set { ProtectedDrawBars = value; }
        }

        public bool DrawDataLabels
        {
            get { return ProtectedDrawDataLabels; }
            set { ProtectedDrawDataLabels = value; }
        }

        public float BarWidth
        {
            get { return ProtectedBarWidth; }
            set { ProtectedBarWidth = value; }
        }

        public float PenWidth
        {
            get { return _penWidth; }
            set { _penWidth = value; }
        }

        public virtual void AddPoint(PointF pt, int index)
        {
            ProtectedPoints.Add(new PointD(pt.X, pt.Y));
            ScreenPoints.Add(SimplePlot.GraphToScreen(pt));
            Indexes.Add(index);
        }

        public virtual void AddPoint(PointD pt, int index)
        {
            ProtectedPoints.Add(pt);
            ScreenPoints.Add(SimplePlot.GraphToScreen(pt));
            Indexes.Add(index);
        }


        public void Clear()
        {
            ProtectedPoints.Clear();
            ScreenPoints.Clear();
            Indexes.Clear();
        }

        public override float GetMinXVal()
        {
            float mn = float.MaxValue;
            foreach (PointD pt in ProtectedPoints)
            {
                if (pt.X < mn)
                    mn = (float) pt.X;
            }
            return mn;
        }

        public override float GetMaxXVal()
        {
            float max = float.MinValue;
            foreach (PointD pt in ProtectedPoints)
            {
                if (pt.X > max)
                    max = (float) pt.X;
            }
            return max;
        }

        public override float GetMinYVal()
        {
            float mn = float.MaxValue;
            foreach (PointD pt in ProtectedPoints)
            {
                if (pt.Y < mn)
                    mn = (float) pt.Y;
            }
            return mn;
        }

        public override float GetMaxYVal()
        {
            float max = 0;
            foreach (PointD pt in ProtectedPoints)
            {
                if (pt.Y > max)
                    max = (float) pt.Y;
            }
            return max;
        }

        public override void Draw(Graphics g)
        {
            float width = ProtectedDotWidth;

            StringFormat format = new StringFormat
                                      {
                                          Alignment = StringAlignment.Center,
                                          LineAlignment = StringAlignment.Center
                                      };

            using (Brush labelBrush = new SolidBrush(Color.Black))
            using (Font labelFont = new Font(FontFamily.GenericSansSerif, 10, FontStyle.Regular))
            using (SolidBrush brush = new SolidBrush(Color))
            using (Pen pen = new Pen(Color, _penWidth))
            {
                if (ProtectedDrawBars)
                {
                    DrawPointsUsingBars(g, brush, format, labelFont, labelBrush);
                }
                else
                {
                    for (int pointIndex = 0; pointIndex < ScreenPoints.Count; pointIndex++)
                    {
                        PointF pt = ScreenPoints[pointIndex];
                        if (width > 1)
                        {
                            g.FillEllipse(brush, pt.X - width/2.0f, pt.Y - width/2.0f, width, width);
                        }

                        if (ProtectedLinkDots && pointIndex > 0)
                        {
                            PointF pt0 = ScreenPoints[pointIndex - 1];
                            g.DrawLine(pen, pt0, pt);
                        }
                    }
                }
            }
        }

        private void DrawPointsUsingBars(Graphics g, SolidBrush brush, StringFormat format, Font labelFont,
                                         Brush labelBrush)
        {
            // If BarWidth <= 0, we'll use whatever bar size makes sense 
            // to connect each consecutive pair of points.  However, we probably don't want to 
            // use very wide bars for points that are unusually far apart.  So, we'll place an upper
            // limit on bar widths:
            float barDxLimit = float.MaxValue;
            if (ProtectedBarWidth <= 0)
            {
                List<float> steps = new List<float>();
                for (int pointIndex = 1; pointIndex < ScreenPoints.Count; pointIndex++)
                {
                    PointF pt = ScreenPoints[pointIndex];
                    float leftDx = 1;
                    PointF ptLeft = ScreenPoints[pointIndex - 1];
                    leftDx = (pt.X - ptLeft.X)/2;
                    steps.Add(leftDx);
                }
                steps.Sort();
                if (steps.Count > 1)
                {
                    barDxLimit = steps[(int) (steps.Count*0.9)];
                }
            }

            for (int pointIndex = 0; pointIndex < ScreenPoints.Count; pointIndex++)
            {
                PointF pt = ScreenPoints[pointIndex];

                // Determine the dX between this point's center and its left and right edges:
                float leftDx = 1;
                float rightDx = 1;
                if (ProtectedBarWidth <= 0)
                {
                    // This bar should extend halfway to the neighboring points...but not farther than BarDxLimit:
                    if (pointIndex > 0)
                    {
                        PointF ptLeft = ScreenPoints[pointIndex - 1];
                        leftDx = (pt.X - ptLeft.X)/2;
                        if (leftDx < 0) leftDx = 1;
                    }
                    if (pointIndex < ScreenPoints.Count - 1)
                    {
                        PointF ptRight = ScreenPoints[pointIndex + 1];
                        rightDx = (ptRight.X - pt.X)/2;
                        if (rightDx < 0) rightDx = 1;
                    }
                    if (pointIndex == 0) leftDx = rightDx;
                    if (pointIndex == ScreenPoints.Count - 1) rightDx = leftDx;
                    if (leftDx > barDxLimit) leftDx = barDxLimit;
                    if (rightDx > barDxLimit) rightDx = barDxLimit;
                }
                else
                {
                    leftDx = SimplePlot.GetPixelsPerX()*ProtectedBarWidth/2;
                    rightDx = leftDx;
                }
                float w = leftDx + rightDx;
                if (w < 1) w = 1;
                if (VerticalGradient)
                {
                    Color startColor = VerticalGradientFrom;
                    Color endColor = VerticalGradientTo;

                    RectangleF rect = new RectangleF(pt.X - leftDx, SimplePlot.ConvertYToPixel(SimplePlot.MaxY), w,
                                                     Math.Abs(SimplePlot.ConvertYToPixel(0) -
                                                              SimplePlot.ConvertYToPixel(SimplePlot.MaxY)) + 1);
                    using (
                        LinearGradientBrush gradBrush = new LinearGradientBrush(rect, startColor, endColor,
                                                                                LinearGradientMode.Vertical))
                    {
                        ColorBlend blend = new ColorBlend
                                               {
                                                   Colors = new[] {startColor, startColor, endColor},
                                                   Positions = new[] {0.0F, 0.74F, 1.0F}
                                               };
                        gradBrush.InterpolationColors = blend;
                        g.FillRectangle(gradBrush, pt.X - leftDx, pt.Y, w,
                                        Math.Abs(SimplePlot.ConvertYToPixel(0) - pt.Y));
                        using (Pen borderPen = new Pen(endColor))
                        {
                            g.DrawRectangle(borderPen, pt.X - leftDx, pt.Y, w,
                                            Math.Abs(SimplePlot.ConvertYToPixel(0) - pt.Y));
                        }
                    }
                }
                else
                {
                    g.FillRectangle(brush, pt.X - leftDx, pt.Y, w, Math.Abs(SimplePlot.ConvertYToPixel(0) - pt.Y));
                }

                if (ProtectedDrawDataLabels)
                {
                    PointD labelPt = Points[pointIndex];
                    float tval = (float) labelPt.Y;
                    string cstr;

                    if (tval == (int) tval)
                    {
                        cstr = tval >= 100000
                                   ? string.Format("{0}K", (int) (0.5 + tval/1000.0))
                                   : string.Format("{0}", (int) tval);
                    }
                    else
                    {
                        if (tval < .01)
                            cstr = string.Format("{0:F3}", tval);
                        else if (tval < 1)
                            cstr = string.Format("{0:F2}", tval);
                        else if (tval < 1000)
                            cstr = string.Format("{0:F1}", tval);
                        else cstr = string.Format("{0:F0}", tval);
                    }
                    SizeF strsize = g.MeasureString(cstr, labelFont);
                    int labelY = (int) pt.Y - (int) strsize.Height - 3;
                    if (labelY < SimplePlot.Padding.TopPadding) labelY = SimplePlot.Padding.TopPadding;
                    Rectangle labelRect = new Rectangle((int) pt.X - (int) strsize.Width/2 - 2,
                                                        labelY,
                                                        (int) strsize.Width + 4, (int) strsize.Height);

                    g.DrawString(cstr, labelFont, labelBrush, labelRect, format);
                }
            }
        }

        public static void DrawRoundedRectangle(Graphics g, Rectangle r, int radius, Pen p)
        {
            GraphicsPath gp =
                new GraphicsPath();

            gp.AddArc(r.X, r.Y, radius, radius, 180, 90);
            gp.AddArc(r.X + r.Width - radius, r.Y, radius, radius, 270, 90);
            gp.AddArc(r.X + r.Width - radius, r.Y + r.Height - radius, radius, radius, 0, 90);
            gp.AddArc(r.X, r.Y + r.Height - radius, radius, radius, 90, 90);
            gp.CloseAllFigures();

            g.DrawPath(p, gp);
        }

        public override void Calculate()
        {
            for (int i = 0; i < ProtectedPoints.Count; i++)
            {
                ScreenPoints[i] = SimplePlot.GraphToScreen(ProtectedPoints[i]);
            }
        }
    }

    public abstract class GraphObject
    {
        #region constructors

        /// <summary>
        ///     The default constructor.  You MUST set the baseGraph property before using the object
        /// </summary>
        protected GraphObject()
        {
        }

        /// <summary>
        ///     Constructor.  Associates this object with the specified base graph,
        ///     and adds this object to that graph control's object list
        /// </summary>
        /// <param name="simplePlot">The parent graph that this object will be associated with</param>
        protected GraphObject(SimplePlot simplePlot)
        {
            SimplePlot = simplePlot;
        }

        #endregion

        #region public properties

        protected SimplePlot SimplePlot = null;

        /// <summary>
        ///     Gets or sets the parent graph that this graph is associated with.  When the setter is called
        ///     this object is removed from the previous base graph's object list, and then added to the new
        ///     base graph's object list.
        /// </summary>
        public SimplePlot PlotObject
        {
            get { return SimplePlot; }
            set
            {
                if (SimplePlot != null) SimplePlot.RemoveObject(this);
                SimplePlot = value;
                SimplePlot.AddObject(this);
            }
        }

        #endregion

        #region abstract methods

        /// <summary>
        ///     Called when the parent graph has been resized, or the axes have been changed.  This object should
        ///     recalculate the pixel locations of any objects it draws on the graph.
        /// </summary>
        public abstract void Calculate();

        /// <summary>
        ///     Called when the parent graph is drawn.  This object sould draw itself in the specified graphics context
        /// </summary>
        /// <param name="g">The graphics context used for drawing the object</param>
        public abstract void Draw(Graphics g);

        #endregion

        #region virtual methods

        /// <summary>
        ///     Called when the parent graph has been asked to autoscale itself.  This is the value of the minimum for
        ///     the x-axis that would display the current graph "optimally".
        /// </summary>
        /// <returns>minimum x-value</returns>
        public virtual float GetMinXVal()
        {
            return 0;
        }

        /// <summary>
        ///     Called when the parent graph has been asked to autoscale itself.  This is the value of the maximum for
        ///     the x-axis that would display the current graph "optimally".
        /// </summary>
        /// <returns></returns>
        public virtual float GetMaxXVal()
        {
            return 0;
        }

        /// <summary>
        ///     Called when the parent graph has been asked to autoscale itself.  This is the value of the minimum for
        ///     the y-axis that would display the current graph "optimally".
        /// </summary>
        /// <returns></returns>
        public virtual float GetMinYVal()
        {
            return 0;
        }

        /// <summary>
        ///     Called when the parent graph has been asked to autoscale itself.  This is the value of the maximum for
        ///     the y-axis that would display the current graph "optimally".
        /// </summary>
        /// <returns></returns>
        public virtual float GetMaxYVal()
        {
            return 0;
        }

        /// <summary>
        ///     Called to save statistics to XML file.
        /// </summary>
        /// <returns></returns>
        public virtual XElement SaveToXElement(string colName, string rowKey)
        {
            return null;
        }

        #endregion
    }

    public class ACGTLegend : GraphObject
    {
        private readonly Color _colorA = Color.Blue;
        private readonly Color _colorC = Color.Blue;
        private readonly Color _colorG = Color.Blue;
        private readonly Color _colorT = Color.Blue;

        private PointF _basePoint;

        public ACGTLegend()
        {
        }

        public ACGTLegend(SimplePlot simplePlot, Color colorA, Color colorC, Color colorG, Color colorT)
            : base(simplePlot)
        {
            SimplePlot = simplePlot;
            SimplePlot.AddObject(this);

            _colorA = colorA;
            _colorC = colorC;
            _colorG = colorG;
            _colorT = colorT;
        }

        public override void Calculate()
        {
            _basePoint = SimplePlot.GraphToScreen(new PointF(SimplePlot.MaxX, SimplePlot.MaxY));
        }

        public override void Draw(Graphics g)
        {
            int xmax = (int) _basePoint.X - 15;
            int xmin = xmax - 40;
            int y = (int) _basePoint.Y;
            const int ystep = 15;

            DrawBaseLegend(g, "A", _colorA, xmin, xmax, y + ystep*1);
            DrawBaseLegend(g, "C", _colorC, xmin, xmax, y + ystep*2);
            DrawBaseLegend(g, "G", _colorG, xmin, xmax, y + ystep*3);
            DrawBaseLegend(g, "T", _colorT, xmin, xmax, y + ystep*4);
        }

        private void DrawBaseLegend(Graphics g, string name, Color color, int xmin, int xmax, int y)
        {
            int width = xmax - xmin + 1;
            int xmid = xmin + width/2;

            Font font = new Font("Microsoft Sans Serif", 9);

            const int dotWidth = 14;

            SizeF size = g.MeasureString(name, font);

            using (SolidBrush brush = new SolidBrush(color))
            using (Pen pen = new Pen(color, 2))
            {
                g.FillEllipse(brush, xmid - dotWidth/2, y - dotWidth/2, dotWidth, dotWidth);
                g.DrawLine(pen, xmin, y, xmax, y);
                g.DrawString(name, font, brush, xmin - size.Width - 3, y - size.Height/2);
            }
        }
    }

    public class AxisLabel
    {
        public float Position;
        public string Value;
    }

    public class BarPlotPoint
    {
        public AxisLabel Label;
        public float XValue;
        public string YLabel;
        public float YValue;
    }

    public class PointD
    {
        public double X;
        public double Y;

        public PointD(double x, double y)
        {
            X = x;
            Y = y;
        }
    }
}