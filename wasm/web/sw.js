// meshseal web app — service worker (PWA offline shell).
//
// Strategy: cache-first for the same-origin static shell, with runtime
// caching of any first-seen same-origin GET. Because the engine
// (meshseal.js + meshseal.wasm) is large and immutable per build, caching
// it makes repeat visits instant and lets the whole app work offline once
// installed.
//
// IMPORTANT: bump CACHE on every deploy (it is the cache-bust key). The new
// service worker takes over on the next load (skipWaiting + clients.claim)
// and purges older caches on activate, so users never get a stale engine.
const CACHE = 'meshseal-v0.1.3-7';

const ASSETS = [
  './',
  './index.html',
  './app.js',
  './worker.js',
  './viewer.js',
  './meshseal.js',
  './meshseal.wasm',
  './manifest.json',
  './icon-192.png',
  './icon-512.png',
  './icon-180.png',
  './favicon-32.png',
  './vendor/three.module.min.js',
  './vendor/STLLoader.js',
  './vendor/OrbitControls.js',
];

self.addEventListener('install', (event) => {
  event.waitUntil(
    caches.open(CACHE)
      .then((cache) => cache.addAll(ASSETS))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.filter((k) => k !== CACHE).map((k) => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', (event) => {
  const req = event.request;
  if (req.method !== 'GET') return;

  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return; // don't touch cross-origin

  event.respondWith(
    caches.match(req, { ignoreSearch: true }).then((hit) => {
      if (hit) return hit;
      return fetch(req).then((resp) => {
        // Runtime-cache successful same-origin responses (covers the very
        // first visit, before install's addAll has populated everything).
        if (resp && resp.ok && resp.type === 'basic') {
          const copy = resp.clone();
          caches.open(CACHE).then((cache) => cache.put(req, copy));
        }
        return resp;
      }).catch(() => {
        // Offline navigation → serve the cached app shell.
        if (req.mode === 'navigate') return caches.match('./index.html');
        return Response.error();
      });
    })
  );
});
