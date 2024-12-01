
# International Standard Name Identifier (ISNI)

This document is a brief analysis of ISNI codes, which are intended to uniquely identify researchers as used in ORCID codes. This analysis is focused entirely on the academic & research community. It will describe various shortcomings of ISNI and will ultimately favour alternative ways to uniquely identify researchers.


1. ISNI, https://isni.org/


## Standards

International Standard Name Identifier (ISNI) is defined by ISO 27729 (2024), which replaces the standard ISO 27729 (2012; Edition 1):
> ISO 27729:2024 Information and documentation —
> International standard name identifier (ISNI)
> https://www.iso.org/standard/87177.html

Unfortunately, the standard is not open source and I did not have any incentive to buy a copy of it (partly for technical reasons that will become more obvious is the next section).


## ISNI vs ISBN

ISNI seems to resemble (at least superficially) the ISBN coding system. Although I did not perform a detailed research on the technical aspects - the available details and the nature of ISNI and of the publishing community which is behind ISNI support this view.

ISBNs are a very simple method to identify uniquely a book. It is not meaningful to have a set of different books under the same ISBN. Each book has assigned a single ISBN code.

Similarly, ISNI identifies one individual (or one organization): unlike ISBNs, it is intended to identify a collection of works published "under"/"by" the same ISNI. The ISBN-model may not be the ideal technical solution for this scenario.

Furthermore, most articles are written by a group of authors: there are rarely only 1 or 2 authors; most often there are more than 4-5 authors. Using 16 digits codes for each author wastes a lot of space. An alternative solution will be presented below.


## Academic Community

### Registrars

Most researchers have their first contact with academic research when they enrol as students in a university. The academic activity increases during the PhD: therefore, either the PhD school and/or the universities should be the primary registrars.

The universities have the ability to check and validate a person as well.

Consider as an example the European Union: each student is either a citizen of a member country (and therefore identifiable in the EU) or needs to apply for a study visa and is therefore also fully identifiable in the Schengen-space. Other registrars do not have any legal framework to identify a person.

It also makes sense that all students across the EU are identifiable as students across the entire Schengen space. As such, an identification system is both necessary and meaningful. Furthermore, it would be easy to extend it to provide unique identification for researchers across the EU.

Technical notes:
- Unique Identifiers: could also encode the region where they were issued, e.g. EU Schengen space, US/Canada, China, etc.
- Other registrars or interested organisations could access specific mechanisms to validate a unique identifier issued in another region; those organisations have most likely access to additional personal information which can be used during the validation steps.
  - 2-Way Identification: the new registrar sends a request to the initial registrar to validate an existing unique identifier;
  - the initial registrar requests a hash computed on a set of personal data (the order of the data and the hashing method specified in the request);
  - the new registrar computes this hash and sends it for validation;
  - the initial registrar validates the unique identifier or rejects the validation;

Researchers and students could receive a robust unique identifier in this way.


## Open Researcher and Contributor ID (ORCID)

The ORCID seems to be based on ISNIs. However, the official suggestion is to use the ISNI as a physical url (with the "https://" part), which hints strongly at a click-bait. The purpose of the unique identifiers is to uniquely identify researchers when they move from one institution to a different one: it is easy to search the literature for the published work, if there is such a unique identifier; the physical link doesn't seem necessary and you do not gain anything new from accessing the ORCID site.

Furthermore, the claim of over 14 million registered users seems also highly unrealistic and indicates a total lack of mechanisms to validate the researchers.

1. Wikipedia: ORCID
> https://en.wikipedia.org/wiki/ORCID
2. GitHub: ORCID code base
> https://github.com/ORCID/


## Alternative Models

**Digital Author Identifier (DAI):** used in the Netherlands.
1. DAI: https://wiki.surfnet.nl/display/standards/DAI
2. Wikipedia: DAI. https://en.wikipedia.org/wiki/Digital_Author_Identifier
3. Author identification: https://repinf.pbworks.com/w/page/13779410/Author%20identification


## Proposal: Academic ID (ACID)

### Character Set

The proposal is to use the upper-case ASCI letters:
- A - Z (26 variants);
- Digits: 1 - 9 could be added in the future; digit 0 is omitted as it could be misidentified with the letter "O".
- Additional characters: can be used inside the string to separate 2 tokens, e.g."-", "+", "\*", "%", "#", "@";
- Validation code: could also use "$" and/or "!";
- Example: "AZYY+AUG_", where "_" stands for the validation code;
- Tokens are formed of 4 characters and are separated by one of the special characters;
- The initial acid codes would be formed by 2 tokens of 4 characters separated by one special character; they could be extended to the left with more tokens.


### Non-Latin Character Sets

Countries using non-Latin character sets could provide mappings to numeric codes or the ASCI set. The characters used should be visually easily identifiable. Simple mappings could be provided for Greek => ASCI, Cyrillic => ASCI, Chinese symbols => ASCI, Japanese Katakana => ASCI, etc. A user having access only to a printed material could easily map manually/visually the non-ASCI codes to their ASCI counterpart.


### Validation Codes

The last character is the validation code. The purpose of the validation code is to ensure the correctness of the entire string:

**Computation:**
- Validation code = sum(position from end \* numeric equivalent of character) (mod p);
- Note: Simple sum: does not identify character inversions; multiplying each character with its position should be more robust.
- Note: gcd(consecutive numbers) = 1;
- Multiplication vector: could be reinitialized for each token with the values {4, 3, 2, 1}; the last token is an exception as it skips the validation character (the last character). Token separators could be multiplied by 5.

**Prime p:**
- Using a prime is more robust: ensures uniqueness even in the multiplicative version;
- p > character set used: the mutation of a single character will always generate an invalid validation code!
- Example: 26 (for A-Z) + 9 (for 1-9) = 35;
- p = 41 is a good candidate, which will accommodate other characters as well.
- The validation code needs therefore 41 characters!

**Encoding:**
- Each character is encoded using a numeric value starting at 1;
- Zero is NOT encoded by any of the characters: this will catch errors due to missing characters.
- Only blank characters are 0: the ACID code can be extended to the left with more characters/tokens, and the old validation codes still remain valid.
- Note: extension to left requires multiplication with position from end of string.
- Token separators: the special characters could be mapped to the same values as some of the "A-Z" characters; or they could be encoded with values beyond "A-Z"; both options will give more variants of IDs.

The simple design (2 tokens of letters; last token only 3 coding characters) offers 26^7 = slightly above 8 billion IDs. Adding 4 characters as separators increases this to 32 billion IDs. Including the digits 1-9 in the character set increases the number to over 257 billion IDs.

A design-variant (of this code) could use the first character to encode the region and/or country where the code was registered.


## TODO

- More analysis...
- Academic ID (ACID): extend/refine new proposal;

