var page = require('webpage').create(),
	system = require('system');

var timeOut = 1000;
if (system.args.length > 3) { timeOut = +system.args[3]; }
page.open(system.args[1], function () {
    window.setTimeout(function () {
        page.render(system.args[2], { quality: '70' });
        phantom.exit();
    }, timeOut);
});
