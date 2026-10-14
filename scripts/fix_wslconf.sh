#!/bin/bash
cat > /etc/wsl.conf <<'EOC'
[boot]
systemd=true
[user]
default=dbs
EOC
cat /etc/wsl.conf
