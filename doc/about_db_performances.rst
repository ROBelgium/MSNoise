.. _aboutdbandperformances:

About Databases and Performances
================================


To quote the SQLite website:

    Appropriate Uses For SQLite

    SQLite is different from most other SQL database engines in that its primary design goal is to be simple:

        * Simple to administer
        * Simple to operate
        * Simple to embed in a larger program
        * Simple to maintain and customize

    Many people like SQLite because it is small and fast. But those qualities are just happy accidents. Users also find that SQLite is very reliable. Reliability is a consequence of simplicity. With less complication, there is less to go wrong. So, yes, SQLite is small, fast, and reliable, but first and foremost, SQLite strives to be simple.

    Simplicity in a database engine can be either a strength or a weakness, depending on what you are trying to do. In order to achieve simplicity, SQLite has had to sacrifice other characteristics that some people find useful, such as high concurrency, fine-grained access control, a rich set of built-in functions, stored procedures, esoteric SQL language features, XML and/or Java extensions, tera- or peta-byte scalability, and so forth. If you need some of these features and do not mind the added complexity that they bring, then SQLite is probably not the database for you. SQLite is not intended to be an enterprise database engine. It is not designed to compete with Oracle or PostgreSQL.

    The basic rule of thumb for when it is appropriate to use SQLite is this: Use SQLite in situations where simplicity of administration, implementation, and maintenance are more important than the countless complex features that enterprise database engines provide. As it turns out, situations where simplicity is the better choice are more common than many people realize.

    Another way to look at SQLite is this: SQLite is not designed to replace Oracle. It is designed to replace fopen(). 
    

To test MSNoise, one can work with a SQLite database. SQLite communication is supported by default in Python (part of the standard library). The major drawback of SQLite is that it doesn't support high concurrency. In the case of MSNoise, this means that only one
Thread (or Process) can interact with the database "at a time". For small batch tests or small runs, that is OK, but when processing larger archives (years of data of 5+ stations), then
the implementation of a MySQL database will allow to process the jobs in parallel.

.. note:: I have been working on some sort of API server layer above a single SQLite database, working as a Queuing system. The API server is the only client of the database, and exchanges data with the code *via* json HTTP requests. Any help, idea, brainstorming on this is welcome!