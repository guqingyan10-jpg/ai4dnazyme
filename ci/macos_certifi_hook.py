"""Provide a bundled CA store for explicitly requested HTTPS input fetching."""
import os
import certifi

os.environ.setdefault("SSL_CERT_FILE", certifi.where())
