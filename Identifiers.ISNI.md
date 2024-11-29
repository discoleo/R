
# International Standard Name Identifier (ISNI)

This document is a brief analysis of ISNI to uniquely identify researchers as used in ORCID codes. This analysis is focused entirely on the academic & research community. It will actually favour alternative ways to uniquely identify researchers.

1. ISNI, https://isni.org/


## Standards

International Standard Name Identifier (ISNI) is defined by ISO 27729 (2024), which replaces the standard ISO 27729 (2012; Edition 1):
> ISO 27729:2024 Information and documentation —
> International standard name identifier (ISNI)
> https://www.iso.org/standard/87177.html

Unfortunately, the standard is not open source and I did not have any incentive to buy a copy of it (partly for technical reasons that will become more obvious is the next section).


## ISNI vs ISBN

ISNI seems to resemble (at least superficially) the ISBN coding system. Although I did not perform a detailed research on the technical aspects - the available details and the nature of ISNI and of the publishing community which is behind ISNI support this view.

ISBNs are a very simple method to identify uniquely a book. It is not meaningful to have a set of different books under the same ISBN.

Similarly, ISNI identifies one individual (or one organization): unlike ISBNs, it is intended to identify a collection of works published "under"/"by" the same ISNI. The ISBN-model may not be the ideal technical solution for this scenario.


## Academic Community

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

> Wikipedia: ORCID
> https://en.wikipedia.org/wiki/ORCID


## TODO

More analysis...
