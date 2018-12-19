# Configuration Examples and Snippets

The filename's prefix defines the configuration's category and shall start with one of the following:

* `section_*.gronf`
  For sections
* `snippet_*.gronf`
  Pure snippets

## Convention

First comment describes the category as, *Section* or *Snippet*, followed by a collon. (the category is optional).

```yaml
%YAML 1.1
# Example: Regional CMT from waveform observations.
--- !grond.Config

# All file paths referenced below are treated relative to the location of this
# configuration file, here we may give a common prefix. E.g. setting it to '..'
# if the configuration file is...
```
