#!/usr/bin/env node
// Incomplete for now

const { createServer } = require('http');

createServer(function (req, res) {
    if (req.url === '/') {
        res.writeHead(200, { 'Content-Type': 'text/html;charset=utf8' });
        res.end(`
        <meta name="color-scheme" content="light dark">
        <div style="width: 800px; margin: 2em auto;">
            <fieldset style="border: 0">
                <legend>Method:</legend>
                <label><input type="radio" id="yallop" checked> Yallop</label><br>
                <label><input type="radio" id="odeh"> Odeh</label>
            </fieldset>
            <fieldset style="border: 0">
                <legend>Time:</legend>
                <label><input type="radio" id="evening" checked> Evening</label><br>
                <label><input type="radio" id="morning"> Morning</label>
            </fieldset>
            <h2>Map</h2>
            Date: <input type="text" value="2022-08-27T00:00:00Z">
            <h2>Calculate</h2>
            <textarea placeholder="2022-08-27T00:00:00Z 21.23 34.2" style="width: 400px; height: 300px"></textarea>
        </div>
        `);
    }
    else res.end('Not supported');
}).listen(4530);
console.log('Open http://127.0.0.1:4530')