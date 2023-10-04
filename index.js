
const fs = require("fs");
// const { WASI } = require("wasi");
// const wasi = new WASI({ version: 'preview1' });
// const importObject = { wasi_snapshot_preview1: wasi.wasiImport };

(async () => {
  const wasm = await WebAssembly.compile(
    fs.readFileSync("index.wasm")
  );
  const wasiShim = { fd_write: console.log, fd_read: console.log, proc_exit: console.log, fd_close: console.log, fd_seek: console.log };
  const instance = await WebAssembly.instantiate(wasm, { wasi_snapshot_preview1: wasiShim });
// console.log(importObject);
//   wasi.start(instance);
})();


// // var utf8Encoder = new TextEncoder("utf8");

// const { WASI } = require("wasi");
// const wasi = new WASI({ version: 'preview1' });
// const importObject = { wasi_snapshot_preview1: wasi.wasiImport };

// WebAssembly.instantiateStreaming(require('fs').promises.readFile('index.wasm'), importObject).then(function (result) {

//   var exports = result.instance.exports;
//   console.log(exports);

// //   fetch('../test/fonts/noto/NotoSans-Regular.ttf').then(function (x) {
// //     return x.arrayBuffer();
// //   }).then(function (fontBlob) {
// //     var heapu8 = new Uint8Array(exports.memory.buffer);
// //     var heapu32 = new Uint32Array(exports.memory.buffer);
// //     var heapi32 = new Int32Array(exports.memory.buffer);

// //     var fontBuffer = exports.malloc(fontBlob.byteLength);
// //     heapu8.set(new Uint8Array(fontBlob), fontBuffer);

// //     var blob = exports.hb_blob_create(fontBuffer, fontBlob.byteLength, 2/*HB_MEMORY_MODE_WRITABLE*/, 0, 0);
// //     var face = exports.hb_face_create(blob, 0);
// //     var font = exports.hb_font_create(face);

// //     var buffer = exports.hb_buffer_create();
// //     {
// //       var text = utf8Encoder.encode('abc');
// //       var text_ptr = exports.malloc(text.byteLength);
// //       heapu8.set(text, text_ptr);
// //       exports.hb_buffer_add_utf8(buffer, text_ptr, text.byteLength, 0, -1);
// //       exports.free(text_ptr);
// //     }
// //     exports.hb_buffer_guess_segment_properties(buffer);
// //     // exports.hb_buffer_set_direction(5); 4: ltr, 5: rtl, 6: ttb. 7: btt

// //     exports.hb_shape(font, buffer, 0, 0);

// //     var length = exports.hb_buffer_get_length(buffer);
// //     var result = [];
// //     var infosPtr32 = exports.hb_buffer_get_glyph_infos(buffer, 0) / 4;
// //     var positionsPtr32 = exports.hb_buffer_get_glyph_positions(buffer, 0) / 4;
// //     var infos = heapu32.subarray(infosPtr32, infosPtr32 + 5 * length);
// //     var positions = heapi32.subarray(positionsPtr32, positionsPtr32 + 5 * length);
// //     for (var i = 0; i < length; ++i) {
// //       result.push({
// //         g: infos[i * 5 + 0],
// //         cl: infos[i * 5 + 2],
// //         ax: positions[i * 5 + 0],
// //         ay: positions[i * 5 + 1],
// //         dx: positions[i * 5 + 2],
// //         dy: positions[i * 5 + 3]
// //       });
// //     }

// //     exports.hb_buffer_destroy(buffer);
// //     exports.hb_font_destroy(font);
// //     exports.hb_face_destroy(face);
// //     exports.hb_blob_create(blob);
// //     exports.free(fontBuffer);

// //     document.body.innerText = JSON.stringify(result);
// //   });
// });