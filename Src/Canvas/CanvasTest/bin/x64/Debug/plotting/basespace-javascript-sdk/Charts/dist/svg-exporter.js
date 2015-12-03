var bs;
(function (bs) {
    bs.svgToPng = function (svgStr, customOptions) {
        var options = {
            cssString: null,
            width: 640,
            height: 480
        };

        var popupWindow, openPopup;
        $.extend(true, options, customOptions);

        openPopup = function () {
            var windowOptions = "width=" + options.width + ",height=" + (options.height + 60) + ",scrollbars=no";
            popupWindow = window.open("", "", windowOptions);
            popupWindow.document.write('<!doctype html><html>');
            popupWindow.document.write("<p>Right click on the image below to download or copy.</p>");
            popupWindow.document.write('<style>#svg-container { position: absolute; top: 0; left: 0; } body { margin: 0; padding: 0; } img { max-width: 95%; } p { font-family: sans-serif; color: #333; padding-top: 5px; margin-bottom: 0; font-size: 13px; text-align: center;}</style>');
            popupWindow.document.write('<img id="png" />');
            popupWindow.document.write('<div id="svg-container" style="width:' + options.width + 'px;height:' + options.height + 'px">' + svgStr + '</div>');
            popupWindow.document.write('<canvas id="canvas-svg" style="width:' + options.width + 'px;height:' + options.height + 'px"></canvas></html>');
            var scriptsLoaded = 0;
            var scriptSources = [
                'https://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js',
                'https://canvg.googlecode.com/svn/trunk/rgbcolor.js',
                'https://canvg.googlecode.com/svn/trunk/canvg.js'];
            var scriptEls;
            var doc = popupWindow.document;

            var onLoad = function () {
                scriptsLoaded++;
                if (scriptsLoaded == scriptSources.length) {
                    var setup = function (cssString, scaleWidth, scaleHeight) {
                        setTimeout(function () {
                            var canvas = $("#canvas-svg");
                            if (cssString)
                                $('#svg-container svg').prepend('<style><![CDATA[' + cssString + ']]></style>');
                            canvg(canvas[0], $("#svg-container").html(), {
                                scaleWidth: scaleWidth,
                                scaleHeight: scaleHeight
                            });
                            var theImage = canvas[0].toDataURL("image/png");

                            $("#png").attr("src", theImage);

                            $('#svg-container, #canvas-svg').remove();
                        }, 250);
                    };
                    script = doc.createElement('script');
                    var cssStringEscaped = options.cssString;
                    cssStringEscaped = cssStringEscaped.replace(/\"/g, "\\\"");
                    script.textContent = "(" + setup.toString() + ")(\"" + cssStringEscaped + "\", '" + options.width + "px', '" + options.height + "px');";
                    doc.body.appendChild(script);
                }
            };

            for (var i = 0; i < scriptSources.length; i++) {
                var script = doc.createElement('script');
                script.src = scriptSources[i];

                script.onload = onLoad;

                doc.head.appendChild(script);
            }
        };

        openPopup();
    };
})(bs || (bs = {}));
