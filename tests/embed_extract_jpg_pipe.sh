#!/usr/bin/env bash

# Copyright 2025
#
# This file is under BSD-3-Clause license.

# Write message using a non-seekable input (FIFO)
echo -e "\nEmbedding a message from FIFO..."

fifo="test-pipe.jpg"
outfile="test-with-message-pipe.jpg"
txtfile="text-jpg-pipe.txt"

cleanup() {
	rm -f "$fifo" "$outfile" "$txtfile"
}
trap cleanup EXIT

rm -f "$fifo"
mkfifo "$fifo" || { echo ERROR; exit 1; }

cat test.jpg > "$fifo" &
cat_pid=$!

../src/outguess -k "secret-key-001" -d message.txt "$fifo" "$outfile" || { echo ERROR; exit 1; }
wait "$cat_pid" || { echo ERROR; exit 1; }

# Retrieve message
echo -e "\nExtracting a message..."
../src/outguess -k "secret-key-001" -r "$outfile" "$txtfile"
cat "$txtfile" | grep "inside of the image" || { echo ERROR; exit 1; }
