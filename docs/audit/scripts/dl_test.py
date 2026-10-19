import ast, os, sys, threading, http.server, socketserver, time, urllib.request, urllib.error, urllib.parse, tempfile, json
src=open('dem-to-stl/sources.py',encoding='utf-8').read()
tree=ast.parse(src)
want={'_http_get','_download'}
ns={'os':os,'time':time,'urllib':urllib,'json':json,'HTTP_TIMEOUT':5,'HTTP_RETRIES':3,'USER_AGENT':'test'}
class SourceUnavailable(RuntimeError): pass
ns['SourceUnavailable']=SourceUnavailable
for node in tree.body:
    if isinstance(node,ast.FunctionDef) and node.name in want:
        exec(compile(ast.Module([node],[]),'sources.py','exec'),ns)
hits={'trunc':0}
class H(http.server.BaseHTTPRequestHandler):
    protocol_version='HTTP/1.1'
    def log_message(self,*a): pass
    def do_GET(self):
        if self.path=='/trunc':
            hits['trunc']+=1
            self.send_response(200); self.send_header('Content-Length','1000'); self.end_headers()
            self.wfile.write(b'II*\x00'+b'x'*496); self.wfile.flush(); self.close_connection=True
        elif self.path=='/xml200':
            body=b'<?xml version="1.0"?><ServiceExceptionReport>serveur surcharge</ServiceExceptionReport>'
            self.send_response(200); self.send_header('Content-Type','text/xml'); self.send_header('Content-Length',str(len(body))); self.end_headers(); self.wfile.write(body)
        else:
            self.send_response(404); self.send_header('Content-Length','0'); self.end_headers()
srv=socketserver.TCPServer(('127.0.0.1',0),H); port=srv.server_address[1]
threading.Thread(target=srv.serve_forever,daemon=True).start()
base=f'http://127.0.0.1:{port}'
print("1) Réponse tronquée (Content-Length annoncé 1000, 500 octets reçus) :")
try:
    ns['_http_get'](base+'/trunc'); print("   -> retour sans erreur")
except Exception as e:
    print(f"   -> exception {type(e).__module__}.{type(e).__name__} (SourceUnavailable ? {isinstance(e,SourceUnavailable)}) ; tentatives du serveur : {hits['trunc']}")
print("2) HTTP 200 avec corps d'erreur XML, via _download (cache) :")
d=tempfile.mkdtemp(); dest=os.path.join(d,'tile.tif')
r1=ns['_download'](base+'/xml200',dest)
print(f"   1er appel : fichier écrit = {os.path.getsize(dest)} octets, début = {open(dest,'rb').read(20)!r}")
srv.shutdown(); r2=ns['_download'](base+'/xml200',dest)
print(f"   2e appel (serveur arrêté) : renvoie le fichier en cache sans le revalider -> {r2 == dest}")
