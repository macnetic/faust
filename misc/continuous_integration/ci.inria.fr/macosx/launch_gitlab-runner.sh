#!/bin/bash

nohup /usr/local/bin/gitlab-runner run --working-directory /Users/ci --config /Users/ci/.gitlab-runner/config.toml --service gitlab-runner --syslog &
